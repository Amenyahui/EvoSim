"""Configuration for current execution and historical result audits."""

from dataclasses import asdict, dataclass

import numpy as np


def species_units(value, sizes, name):
    values = np.asarray(value, dtype=float)
    if values.ndim == 0:
        values = np.full(len(sizes), float(values))
    if values.shape != (len(sizes),) or not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError(f"{name} must be positive finite scalar or one value per species")
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
            "resource and demand must be positive finite numeric values, not strings/bools"
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
    sP0: float
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
            "sP0",
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
        if self.mapping not in ("identity", "mixed"):
            raise ValueError("mapping must be identity or mixed; use B for a custom sparse map")
        if type(self.shuffle_voc) is not bool:
            raise ValueError("shuffle_voc must be boolean")
        B = np.asarray(self.B, float)
        if (
            B.ndim != 2
            or B.shape[1] != self.n
            or B.shape[0] < 1
            or not np.isfinite(B).all()
            or (B < 0).any()
        ):
            raise ValueError("B must be a finite nonnegative L by K matrix")
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
        for name in ("rho_voc",):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        if self.voc_cost < 0:
            raise ValueError("sigma_Z and voc_cost must be nonnegative")
        normalized_bounds(self.voc_bounds, "voc_bounds")
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
        if self.insect_life_history != "self_feeding_pre_reproductive":
            raise ValueError("explicit self_feeding_pre_reproductive life history required")
        if type(self.sI0) not in (int, float) or not np.isfinite(self.sI0) or self.sI0 != 0:
            raise ValueError("self-feeding insects require sI0=0; survival equals feeding gain")
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
    plant_density_stage="offspring_recruitment",
    insect_independent_density_factor=False,
    immigration=False,
    seed_bank=False,
    overlapping_generations=False,
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
    model_version: object
    life_history: object
    plant_offspring_lambda: object
    insect_offspring_lambda: object
    plant_competition_alpha: object
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
            self.model_version not in {"0.1.2", "0.1.3"}
            or self.life_history != "nonoverlapping_sexual_generations"
        ):
            raise ValueError("explicit supported model version and sexual generations required")
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
        provenance = self.parameter_provenance
        required = {
            "plant_offspring_lambda": ps,
            "insect_offspring_lambda": hs,
            "plant_competition_alpha": ps,
        }
        if not isinstance(provenance, dict) or not set(required) | {
            "inherited_ecology",
            "resources",
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
                if field == "plant_competition_alpha":
                    if not isinstance(group[key], dict) or set(group[key]) != set(ps):
                        raise ValueError("alpha provenance needs every recipient/competitor pair")
                    entries.extend(group[key].values())
                else:
                    entries.append(group[key])
        entries += [provenance["inherited_ecology"], provenance["resources"]]
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
        if self.record_observations:
            panel = self.observation_reference_traits
            shapes = {"D": (len(ps), self.n), "plant_voc_independent_trait": (len(ps), len(B))}
            if not isinstance(panel, dict) or set(panel) != set(shapes):
                raise ValueError("observation_reference_traits requires a fixed D/Z panel")
            for name, shape in shapes.items():
                values = np.asarray(panel[name])
                if (
                    values.dtype.kind not in "iuf"
                    or values.shape != shape
                    or not np.isfinite(values).all()
                    or ((values < 0) | (values > 1)).any()
                ):
                    raise ValueError("invalid observation reference panel")
        from .demography import validate_sexual
        from .spatial import validate_spatial

        validate_spatial(self)
        validate_sexual(self)
        from .genetics import validate

        validate(self, len(B))
        return B, w

    def resolved_demography(self):
        ps, hs = names("P", self.plant_sizes), names("H", self.insect_sizes)
        lp = species_map(self.plant_offspring_lambda, ps, "plant_offspring_lambda", True)
        li = species_map(self.insect_offspring_lambda, hs, "insect_offspring_lambda", True)
        alpha = self.plant_competition_alpha
        if not isinstance(alpha, dict) or set(alpha) != set(ps):
            raise ValueError("plant_competition_alpha needs a complete species matrix")
        return lp, li, np.stack([species_map(alpha[s], ps, "alpha row") for s in ps])


def config_payload(c):
    """Serialize the sole current model contract."""
    return asdict(c)
