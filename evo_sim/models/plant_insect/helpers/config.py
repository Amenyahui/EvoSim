"""Public v1.3 configuration and validation.

The private validation layers preserve established boundary behavior. Only
Config is part of the public API; all scientific modules consume that object.
"""

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


@dataclass
class _TraitConfig:
    plant_sizes: tuple = (80, 80)
    insect_pop: int = 160
    n: int = 2
    generations: int = 200
    record_every: int = 10
    seed: int = 20261019
    mutation_rate: float = 1e-5
    mutation_sd: float = 0.2
    mu_A: float = 1e-5
    sigma_A: float = 0.2
    plant_cost: float = 0.125
    insect_cost: float = 0.125
    sP0: float = 0.1
    sI0: float = 0.1
    candidate_size: int = 8
    olfaction_mode: str = "template"
    shuffle_voc: bool = False
    mapping: str = "identity"
    B: object = None
    weights: object = None
    sigma: float = 0.2
    b: float = 1.0
    tau: float = 1.0
    outside_utility: float = 0.0
    unmatched_insect_fitness: float = 0.0
    # Optional explicit species centers; None retains Speed's U(0,.01) start.
    plant_centers: object = None
    initial_spread: float = 0.01
    insect_initial_max: float = 0.01

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
            "mutation_rate",
            "mu_A",
            "sP0",
            "sI0",
            "unmatched_insect_fitness",
            "initial_spread",
            "insect_initial_max",
        ):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        for name in ("mutation_sd", "sigma_A", "plant_cost", "insect_cost", "tau"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be nonnegative")
        if self.sigma <= 0 or self.olfaction_mode not in ("template", "neutral"):
            raise ValueError("positive sigma and template/neutral mode required")
        if self.mapping not in ("identity", "mixed"):
            raise ValueError("mapping must be identity or mixed; use B for a custom sparse map")
        B = (
            np.asarray(self.B, float)
            if self.B is not None
            else (
                np.eye(self.n)
                if self.mapping == "identity"
                else 0.5 * np.eye(self.n) + 0.5 * np.ones((self.n, self.n)) / self.n
            )
        )
        if (
            B.ndim != 2
            or B.shape[1] != self.n
            or B.shape[0] < 1
            or not np.isfinite(B).all()
            or (B < 0).any()
        ):
            raise ValueError("B must be a finite nonnegative L by K matrix")
        w = (
            np.ones(B.shape[0]) / B.shape[0]
            if self.weights is None
            else np.asarray(self.weights, float)
        )
        if (
            w.shape != (B.shape[0],)
            or not np.isfinite(w).all()
            or (w < 0).any()
            or not np.isclose(w.sum(), 1)
        ):
            raise ValueError("weights must be nonnegative, length L and sum to 1")
        if self.plant_centers is not None:
            centers = np.asarray(self.plant_centers, float)
            if (
                centers.shape != (len(self.plant_sizes), self.n)
                or not np.isfinite(centers).all()
                or ((centers < 0) | (centers > 1)).any()
            ):
                raise ValueError("plant_centers must be species by K, in [0,1]")
        return B, w / w.sum()


@dataclass
class _VOCConfig(_TraitConfig):
    rho_voc: float = 1.0
    mu_Z: float = 1e-5
    sigma_Z: float = 0.2
    voc_cost: float = 0.0
    voc_init_mode: str = "matched"
    choice_mode: str = "outside"

    def validate(self):
        B, w = super().validate()
        for name in ("rho_voc", "mu_Z"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        if self.sigma_Z < 0 or self.voc_cost < 0:
            raise ValueError("sigma_Z and voc_cost must be nonnegative")
        if self.voc_init_mode not in ("matched", "independent"):
            raise ValueError("voc_init_mode must be matched or independent")
        if self.choice_mode not in ("outside", "forced"):
            raise ValueError("choice_mode must be outside or forced")
        return B, w


@dataclass
class _CommunityConfig(_VOCConfig):
    insect_pop: int | None = None
    insect_sizes: tuple = (80, 80)
    resource_capacity: object = 1.0
    max_feeding_demand: object = 1.0
    insect_I_initial_max: object = None
    insect_A_initial_bounds: object = None

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
            raise ValueError("v1.2 retains zero fitness for unmatched insects")
        maximum = (
            np.full(len(self.insect_sizes), self.insect_initial_max)
            if self.insect_I_initial_max is None
            else np.asarray(self.insect_I_initial_max, float)
        )
        if (
            maximum.shape != (len(self.insect_sizes),)
            or not np.isfinite(maximum).all()
            or ((maximum < 0) | (maximum > 1)).any()
        ):
            raise ValueError("insect_I_initial_max must give one [0,1] maximum per insect species")
        bounds = (
            np.tile([0.0, 1.0], (len(self.insect_sizes), 1))
            if self.insect_A_initial_bounds is None
            else np.asarray(self.insect_A_initial_bounds, float)
        )
        if (
            bounds.shape != (len(self.insect_sizes), 2)
            or not np.isfinite(bounds).all()
            or ((bounds < 0) | (bounds > 1)).any()
            or (bounds[:, 0] > bounds[:, 1]).any()
        ):
            raise ValueError(
                "insect_A_initial_bounds must be one valid [lower,upper] pair per species"
            )
        return B, w


@dataclass
class _ResourceConfig(_CommunityConfig):
    # Explicit None rejects un-migrated old configs instead of silently adopting beta=.1.
    sI0: object = None
    insect_life_history: object = None
    resource_capacity: object = None
    max_feeding_demand: object = None
    resource_unit: object = None
    demand_unit: object = None
    feeding_window: object = None
    demand_window: object = None

    def validate(self):
        if self.insect_life_history != "obligate_host_larval_feeding":
            raise ValueError("explicit obligate_host_larval_feeding life history required")
        if type(self.sI0) not in (int, float) or not np.isfinite(self.sI0) or self.sI0 != 0:
            raise ValueError(
                "v1.2.1 requires explicit sI0=0 (beta=0); old .1 configs must not be reused"
            )
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
        # ponytail: reuse the complete historical validator on a resolved copy; no alias state.
        old = _CommunityConfig(
            **{k: v for k, v in asdict(self).items() if k in _CommunityConfig.__dataclass_fields__}
        )
        old.resource_capacity = R
        old.max_feeding_demand = dmax
        result = old.validate()
        self.insect_pop = old.insect_pop
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


def species_map(value, keys, name, positive=False):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise ValueError(f"{name} requires exactly {keys}")
    return np.array([numeric(value[k], name, positive) for k in keys], float)


@dataclass
class Config(_ResourceConfig):
    model_version: object = None
    life_history: object = None
    plant_offspring_lambda: object = None
    insect_offspring_lambda: object = None
    plant_competition_alpha: object = None
    demography: object = None
    illustrative_only: object = None
    parameter_provenance: object = None
    study_site: object = None
    max_individuals_safety: object = None
    evolve_plants: bool = True
    evolve_insects: bool = True
    founder_traits: object = None
    record_observations: bool = True
    spatial: object = None
    genetics: object = None

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
            self.model_version not in ("1.3", "1.4", "1.5", "1.6")
            or self.life_history != "nonoverlapping_generations_obligate_host_larva"
        ):
            raise ValueError("explicit v1.3 and nonoverlapping obligate-host life history required")
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
        for name in ("evolve_plants", "evolve_insects", "record_observations"):
            if type(getattr(self, name)) is not bool:
                raise ValueError(f"{name} must be a boolean")
        if self.model_version == "1.3" and (
            not self.evolve_plants or not self.evolve_insects or self.founder_traits is not None
        ):
            raise ValueError("Evolution interventions require model_version=1.4")
        if not (self.evolve_plants and self.evolve_insects) and self.founder_traits is None:
            raise ValueError("Frozen controls require explicit homogeneous founder_traits")
        if self.founder_traits is not None:
            shapes = {
                "D": (len(ps), self.n),
                "I": (len(hs), self.n),
                "A": (len(hs), len(B)),
                "plant_voc_independent_trait": (len(ps), len(B)),
            }
            if not isinstance(self.founder_traits, dict) or set(self.founder_traits) != set(shapes):
                raise ValueError("founder_traits requires D, I, A, plant_voc_independent_trait")
            for name, shape in shapes.items():
                values = np.asarray(self.founder_traits[name])
                if (
                    values.dtype.kind not in "iuf"
                    or values.shape != shape
                    or not np.isfinite(values).all()
                    or ((values < 0) | (values > 1)).any()
                ):
                    raise ValueError(
                        f"founder_traits {name} requires finite [0,1] values of shape {shape}"
                    )
        if (
            self.model_version in ("1.4", "1.5", "1.6")
            and self.record_observations
            and self.founder_traits is None
        ):
            raise ValueError(
                "Standardized observations require explicit founder_traits reference panel"
            )
        if self.model_version in ("1.5", "1.6"):
            from .spatial import validate_spatial

            validate_spatial(self)
        elif self.spatial is not None:
            raise ValueError("Spatial configuration requires model_version=1.5 or 1.6")
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
    """Preserve the v1.3 serialization schema for unmodified legacy runs."""
    values = asdict(c)
    if c.genetics is None or c.genetics == {"mode": "direct_trait_legacy"}:
        values.pop("genetics")
    if c.model_version not in ("1.5", "1.6"):
        values.pop("spatial")
    if c.model_version == "1.3":
        for name in ("evolve_plants", "evolve_insects", "founder_traits", "record_observations"):
            values.pop(name)
    return values
