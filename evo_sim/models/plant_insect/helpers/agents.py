"""Individual agents stored as aligned arrays, one row per individual.

Species labels annotate individuals; they are not population averages. Arrays
are owned by the current generation. Snapshots copy them before reproduction.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class PlantPopulation:
    ids: np.ndarray
    species: np.ndarray
    defense: np.ndarray
    voc_trait: np.ndarray | None = None
    cells: np.ndarray | None = None
    genotypes: np.ndarray | None = None

    chemical_traits: np.ndarray | None = None

    def __len__(self):
        return len(self.ids)


@dataclass(slots=True)
class InsectPopulation:
    ids: np.ndarray
    species: np.ndarray
    resistance: np.ndarray
    preference: np.ndarray
    cells: np.ndarray | None = None
    genotypes: np.ndarray | None = None

    sex: np.ndarray | None = None  # 0=female, 1=male; assigned once at birth
    birth_host_id: np.ndarray | None = None
    birth_host_species: np.ndarray | None = None  # -1 = unknown founder history
    adult_cohort: np.ndarray | None = None

    def __len__(self):
        return len(self.ids)


@dataclass(slots=True)
class CommunityState:
    plants: PlantPopulation
    insects: InsectPopulation

    def validate(self):
        for pop in (self.plants, self.insects):
            if (
                len(np.unique(pop.ids)) != len(pop)
                or pop.genotypes.shape[:2] != (len(pop), 2)
                or pop.cells.shape != (len(pop), 2)
                or len(pop.species) != len(pop)
            ):
                raise ValueError("unaligned population state")
        if (
            self.insects.sex.shape != (len(self.insects),)
            or not np.isin(self.insects.sex, [0, 1]).all()
        ):
            raise ValueError("insect sex must be aligned female=0/male=1")
        for name, values in self.arrays().items():
            n = (
                len(self.plants)
                if name in ("D", "plant_voc_independent_trait") or name.startswith("plant_")
                else len(self.insects)
            )
            if len(values) != n:
                raise ValueError(f"unaligned {name}")

    def arrays(self):
        """Stable v1.3 on-disk names; no copying or random draws."""
        values = dict(
            D=self.plants.defense,
            I=self.insects.resistance,
            A=self.insects.preference,
            plant_species=self.plants.species,
            insect_species=self.insects.species,
            plant_ids=self.plants.ids,
            insect_ids=self.insects.ids,
            insect_sex=self.insects.sex,
        )
        for name in ("birth_host_id", "birth_host_species", "adult_cohort"):
            value = getattr(self.insects, name)
            if value is not None:
                values["insect_" + name] = value
        if self.plants.chemical_traits is not None:
            values["plant_chemical_traits"] = self.plants.chemical_traits
        else:
            values["plant_voc_independent_trait"] = self.plants.voc_trait
        if self.plants.cells is not None:
            values["plant_cells"] = self.plants.cells
            values["insect_cells"] = self.insects.cells
        if self.plants.genotypes is not None:
            values["plant_genotypes"] = self.plants.genotypes
            values["insect_genotypes"] = self.insects.genotypes
        return values

    def snapshot(self, generation):
        return dict(
            generation=np.array(generation),
            **{key: value.copy() for key, value in self.arrays().items()},
        )

    @classmethod
    def from_arrays(cls, values):
        """Construct from internally validated initial or offspring arrays."""
        return cls(
            PlantPopulation(
                values["plant_ids"],
                values["plant_species"],
                values["D"],
                values.get("plant_voc_independent_trait"),
                values.get("plant_cells"),
                values.get("plant_genotypes"),
                values.get("plant_chemical_traits"),
            ),
            InsectPopulation(
                values["insect_ids"],
                values["insect_species"],
                values["I"],
                values["A"],
                values.get("insect_cells"),
                values.get("insect_genotypes"),
                values["insect_sex"],
                values.get("insect_birth_host_id"),
                values.get("insect_birth_host_species"),
                values.get("insect_adult_cohort"),
            ),
        )
