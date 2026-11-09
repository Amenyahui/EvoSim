"""Configuration for the current model protocol."""

from dataclasses import asdict, dataclass

import numpy as np


def species_units(value, sizes, name, allow_zero=False):
    values = np.asarray(value, dtype=float)
    if values.ndim == 0:
        values = np.full(len(sizes), float(values))
    if (
        values.shape != (len(sizes),)
        or not np.isfinite(values).all()
        or (values < 0 if allow_zero else values <= 0).any()
    ):
        raise ValueError(f"{name} must be finite and nonnegative (demand strictly positive)")
    expanded = np.repeat(values, sizes)
    with np.errstate(over="ignore"):
        total = expanded.sum()
    if not np.isfinite(total):
        raise ValueError(f"{name} total exceeds floating-point range")
    return expanded


def resource_values(value, sizes, prefix):
    if isinstance(value, dict):
        names = [f"{prefix}{s + 1}" for s in range(len(sizes))]
        if set(value) != set(names):
            raise ValueError(f"{prefix} resource mapping must contain exactly {names}")
        value = [value[name] for name in names]
    values = np.asarray(value)
    if values.dtype.kind not in "iuf" or not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError(
            "resource must be nonnegative and demand positive finite numeric values, not strings/bools"
        )
    species_units(value, sizes, prefix)
    return value


@dataclass(kw_only=True)
class _TraitConfig:
    # All values come from YAML (or an explicitly supplied Config); no scientific defaults.
    plant_sizes: tuple
    insect_pop: int | None
    n: int
    generations: int
    record_every: int
    seed: int
    plant_cost: float
    insect_cost: float
    sI0: float
    candidate_size: int
    olfaction_mode: str
    shuffle_voc: bool
    mapping: str
    B: object
    weights: object
    sigma: float
    b: float
    tau: float
    outside_utility: float
    unmatched_insect_fitness: float

    def validate(self):
        for name in ("insect_pop", "n", "generations", "record_every", "candidate_size"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if not self.plant_sizes or any(type(x) is not int or x < 1 for x in self.plant_sizes):
            raise ValueError("plant_sizes must contain positive integers")
        if type(self.seed) is not int or not 0 <= self.seed < 2**32:
            raise ValueError("seed must be an integer in [0,2**32)")
        for name, value in asdict(self).items():
            if isinstance(value, float) and not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
        for name in (
            "sI0",
            "unmatched_insect_fitness",
        ):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        for name in ("plant_cost", "insect_cost", "tau"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be nonnegative")
        if self.sigma <= 0 or self.olfaction_mode not in ("template", "neutral"):
            raise ValueError("positive sigma and template/neutral mode required")
        if type(self.shuffle_voc) is not bool:
            raise ValueError("shuffle_voc must be boolean")
        if self.mapping != "independent_psm_voc" or self.n != 2 or self.B is not None:
            raise ValueError("independent PSM/VOC modules require n=2 and B=null")
        B = np.empty((2, 0))  # Shape carrier; the current mapping has no defense-to-odor matrix.
        w = np.asarray(self.weights, float)
        if (
            w.shape != (B.shape[0],)
            or not np.isfinite(w).all()
            or (w < 0).any()
            or not np.isclose(w.sum(), 1)
        ):
            raise ValueError("weights must be nonnegative, length L and sum to 1")
        return B, w / w.sum()


@dataclass(kw_only=True)
class _VOCConfig(_TraitConfig):
    rho_voc: float
    voc_cost: float
    voc_bounds: object
    choice_mode: str

    def validate(self):
        B, w = super().validate()
        if self.rho_voc is not None or self.voc_bounds is not None:
            raise ValueError("shared modules retire rho_voc/voc_bounds; set them to null")
        if self.voc_cost < 0:
            raise ValueError("voc_cost must be nonnegative")
        if self.choice_mode not in ("outside", "forced"):
            raise ValueError("choice_mode must be outside or forced")
        return B, w


@dataclass(kw_only=True)
class _CommunityConfig(_VOCConfig):
    insect_sizes: tuple
    resource_capacity: object
    max_feeding_demand: object

    def validate(self):
        if not self.insect_sizes or any(type(x) is not int or x < 1 for x in self.insect_sizes):
            raise ValueError("insect_sizes must contain positive integer species counts")
        total = sum(self.insect_sizes)
        if self.insect_pop is not None and (
            type(self.insect_pop) is not int or self.insect_pop != total
        ):
            raise ValueError("insect_pop must equal sum(insect_sizes); omit it to derive the total")
        self.insect_pop = total
        B, w = super().validate()
        species_units(self.resource_capacity, self.plant_sizes, "resource_capacity")
        species_units(self.max_feeding_demand, self.insect_sizes, "max_feeding_demand")
        if self.unmatched_insect_fitness != 0:
            raise ValueError("unmatched insects must have zero ecological performance")
        return B, w


@dataclass(kw_only=True)
class _ResourceConfig(_CommunityConfig):
    insect_life_history: object
    resource_unit: object
    demand_unit: object
    feeding_window: object
    demand_window: object

    def validate(self):
        if self.insect_life_history != "adult_oviposition_larval_feeding":
            raise ValueError("current life history requires adult oviposition and larval feeding")
        if type(self.sI0) not in (int, float) or not np.isfinite(self.sI0) or self.sI0 != 0:
            raise ValueError("larval survival requires sI0=0; no additive survival floor")
        units = ("normalized_edible_dry_biomass", "g_edible_dry_biomass", "mg_edible_dry_biomass")
        if self.resource_unit not in units or self.demand_unit != self.resource_unit:
            raise ValueError("resource/demand require the same recognized edible dry mass unit")
        if (
            not isinstance(self.feeding_window, str)
            or not self.feeding_window.strip()
            or self.demand_window != self.feeding_window
        ):
            raise ValueError("resource/demand require the same explicit feeding window")
        R = resource_values(self.resource_capacity, self.plant_sizes, "P")
        dmax = resource_values(self.max_feeding_demand, self.insect_sizes, "H")
        resolved = _CommunityConfig(
            **{k: v for k, v in asdict(self).items() if k in _CommunityConfig.__dataclass_fields__}
        )
        resolved.resource_capacity = R
        resolved.max_feeding_demand = dmax
        result = resolved.validate()
        self.insect_pop = resolved.insect_pop
        return result


DEMOGRAPHY = dict(
    plant_survival="bernoulli",
    insect_survival="bernoulli",
    offspring_distribution_given_survival="poisson",
    plant_density_stage="vacant_microsite_lottery",
    insect_independent_density_factor=False,
    immigration=False,
    seed_bank=False,
    plant_overlapping_generations=True,
    insect_overlapping_generations=False,
)


def names(prefix, sizes):
    return [f"{prefix}{s + 1}" for s in range(len(sizes))]


def numeric(value, name, positive=False):
    if (
        type(value) not in (int, float)
        or not np.isfinite(value)
        or (value <= 0 if positive else value < 0)
    ):
        raise ValueError(
            f"{name} must be finite numeric and " + ("positive" if positive else "nonnegative")
        )
    return value


def normalized_bounds(value, name):
    bounds = np.asarray(value)
    if (
        bounds.dtype.kind not in "iuf"
        or bounds.shape != (2,)
        or not np.isfinite(bounds).all()
        or not 0 <= bounds[0] < bounds[1] <= 1
    ):
        raise ValueError(f"{name} requires [lower, upper] within normalized [0,1]")
    return bounds


def species_map(value, keys, name, positive=False):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise ValueError(f"{name} requires exactly {keys}")
    return np.array([numeric(value[k], name, positive) for k in keys], float)


@dataclass(kw_only=True)
class Config(_ResourceConfig):
    functional_modules: object = None
    chemical_processing: object
    contact_response: object
    model_version: object
    life_history: object
    plant_offspring_lambda: object
    insect_offspring_lambda: object
    plant_demography: object
    demography: object
    illustrative_only: object
    parameter_provenance: object
    study_site: object
    max_individuals_safety: object
    record_observations: bool
    spatial: object
    genetics: object
    sexual_reproduction: object
    observation_reference_traits: object

    def validate(self):
        for sizes in (self.plant_sizes, self.insect_sizes):
            if not sizes or any(type(x) is not int or x < 0 for x in sizes):
                raise ValueError(
                    "initial sizes must contain nonnegative integer counts; keep every species slot"
                )
        if self.insect_pop is not None and (
            type(self.insect_pop) is not int or self.insect_pop != sum(self.insect_sizes)
        ):
            raise ValueError("insect_pop, if supplied, must equal initial sum(insect_sizes)")
        if type(self.generations) is not int or self.generations < 0 or self.record_every != 1:
            raise ValueError(
                "nonnegative generations and record_every=1 required for full demographic audit"
            )
        if (
            self.model_version != "0.1.12"
            or self.life_history != "resource_coupled_plants_larval_insects"
            or not isinstance(self.genetics, dict)
            or self.genetics.get("schema_version") != "slim-qtl-6"
        ):
            raise ValueError(
                "current model_version=0.1.12 and slim-qtl-6 resource-coupled plant life history required"
            )
        if self.demography != DEMOGRAPHY or any(
            type(self.demography[k]) is not type(v) for k, v in DEMOGRAPHY.items()
        ):
            raise ValueError(
                "explicit Bernoulli/Poisson demography flags required; no extra density/rescue processes"
            )
        if (
            type(self.max_individuals_safety) is not int
            or not 1 <= self.max_individuals_safety < 2**63
        ):
            raise ValueError("explicit positive int64 max_individuals_safety required")
        if (
            type(self.illustrative_only) is not bool
            or not isinstance(self.study_site, str)
            or not self.study_site.strip()
        ):
            raise ValueError("explicit illustrative_only and study_site required")
        ps, hs = names("P", self.plant_sizes), names("H", self.insect_sizes)
        self.resolved_demography()
        self.validate_plant_demography()
        self.validate_ecology()
        provenance = self.parameter_provenance
        required = {
            "plant_offspring_lambda": ps,
            "insect_offspring_lambda": hs,
        }
        if not isinstance(provenance, dict) or not set(required) | {
            "inherited_ecology",
            "resources",
            "plant_demography",
        } <= set(provenance):
            raise ValueError(
                "parameter_provenance must describe every lambda/alpha and inherited ecology/resources"
            )
        entries = []
        for field, keys in required.items():
            group = provenance[field]
            if not isinstance(group, dict) or set(group) != set(keys):
                raise ValueError(f"provenance {field} species keys mismatch")
            for key in keys:
                entries.append(group[key])
        entries += [
            provenance[key] for key in ("inherited_ecology", "resources", "plant_demography")
        ]
        for entry in entries:
            if (
                not isinstance(entry, dict)
                or entry.get("kind") not in ("measured", "sourced", "theoretical_placeholder")
                or not isinstance(entry.get("source"), str)
                or not entry["source"].strip()
            ):
                raise ValueError("each parameter provenance requires kind and nonempty source")
            if entry["kind"] == "theoretical_placeholder" and not self.illustrative_only:
                raise ValueError("theoretical placeholders require illustrative_only=true")
        # Validate unchanged ecology through a small positive-count proxy; do not alter real zero species.
        args = {k: v for k, v in asdict(self).items() if k in _ResourceConfig.__dataclass_fields__}
        args.update(
            plant_sizes=[1] * len(ps),
            insect_sizes=[1] * len(hs),
            insect_pop=len(hs),
            generations=max(1, self.generations),
        )
        B, w = _ResourceConfig(**args).validate()
        if type(self.record_observations) is not bool:
            raise ValueError("record_observations must be boolean")
        from evo_sim.models.plant_insect.biology.traits import validate as validate_modules

        validate_modules(self.functional_modules)
        if self.record_observations:
            panel = self.observation_reference_traits
            if not isinstance(panel, dict) or set(panel) != {
                "plant_chemical_traits",
                "reference_ids",
            }:
                raise ValueError("reference panel requires chemical traits and reference_ids")
            ids = panel["reference_ids"]
            if (
                not isinstance(ids, list)
                or not ids
                or any(not isinstance(x, str) or not x.strip() for x in ids)
                or len(set(ids)) != len(ids)
            ):
                raise ValueError("unique nonempty reference IDs required")
            values = np.asarray(panel["plant_chemical_traits"])
            if (
                values.dtype.kind not in "iuf"
                or values.shape != (len(ids), 5)
                or not np.isfinite(values).all()
                or ((values < 0) | (values > 1)).any()
            ):
                raise ValueError("invalid observation reference panel")
        from evo_sim.models.plant_insect.biology.insects.demography import validate_sexual
        from evo_sim.models.plant_insect.biology.spatial import validate_spatial

        validate_spatial(self)
        validate_sexual(self)
        from evo_sim.models.plant_insect.biology.genetics.qtl import validate

        validate(self, len(B))
        return B, w

    def resolved_demography(self):
        ps, hs = names("P", self.plant_sizes), names("H", self.insect_sizes)
        lp = species_map(self.plant_offspring_lambda, ps, "plant_offspring_lambda", True)
        li = species_map(self.insect_offspring_lambda, hs, "insect_offspring_lambda", True)
        return lp, li

    def validate_plant_demography(self):
        spec = self.plant_demography
        rates = (
            "background_mortality_rate",
            "defoliation_mortality_rate",
            "resource_recovery_rate",
            "tolerance_recovery_rate",
        )
        required = set(rates) | {
            "newborn_resource_fraction",
            "microsites_per_cell",
        }
        if not isinstance(spec, dict) or set(spec) != required:
            raise ValueError(f"plant_demography requires exactly {sorted(required)}")
        for key in rates:
            numeric(spec[key], key)
        for key in ("microsites_per_cell",):
            if type(spec[key]) is not int or not 1 <= spec[key] < 2**63:
                raise ValueError(f"{key} must be a positive int64 integer")
        numeric(spec["newborn_resource_fraction"], "newborn_resource_fraction")
        if spec["newborn_resource_fraction"] > 1:
            raise ValueError("newborn_resource_fraction must be in [0,1]")

    def validate_ecology(self):
        """Require explicit theoretical dose, nutritional and contact parameters."""
        processing = self.chemical_processing
        vectors = {"dose_scale", "processing_scale", "toxicity_weights", "antinutrition_weights"}
        if not isinstance(processing, dict) or set(processing) != vectors | {
            "nutrition_exponent",
            "provenance",
        }:
            raise ValueError("complete chemical_processing parameters required")
        for name in vectors:
            values = np.asarray(processing[name])
            if (
                values.shape != (2,)
                or values.dtype.kind not in "iuf"
                or not np.isfinite(values).all()
                or (values < 0).any()
            ):
                raise ValueError(
                    f"chemical_processing.{name} requires two nonnegative finite coefficients"
                )
        numeric(processing["nutrition_exponent"], "nutrition_exponent")
        contact = self.contact_response
        if not isinstance(contact, dict) or set(contact) != {"intercept", "slope", "provenance"}:
            raise ValueError("contact_response requires intercept, slope and provenance")
        if type(contact["intercept"]) not in (int, float) or not np.isfinite(contact["intercept"]):
            raise ValueError("finite contact intercept required")
        numeric(contact["slope"], "contact slope")
        for spec in (processing, contact):
            provenance = spec["provenance"]
            if (
                not isinstance(provenance, dict)
                or provenance.get("kind") != "theoretical_placeholder"
                or not isinstance(provenance.get("source"), str)
                or not provenance["source"].strip()
                or not self.illustrative_only
            ):
                raise ValueError("new ecological equations require explicit theoretical provenance")


def config_payload(c):
    """Serialize the sole current model contract."""
    return asdict(c)
