# EvoSim — Project Instructions for Codex

> 本文件是 EvoSim 项目的项目级开发规范。所有 Codex 开发任务都必须遵守本文件。
>
> 如果具体开发任务与本文件冲突，必须先指出冲突并征求确认，不得自行决定修改科学模型。

---

# 1. Project Overview

## 1.1 Project Identity

- **Project:** EvoSim
- **Repository:** https://github.com/Amenyahui/EvoSim
- **Research Domain:** Evolutionary Biology, Theoretical Ecology, Quantitative Genetics, Agent-Based Modeling.

EvoSim is a theoretical, individual-based evolutionary simulation framework designed to investigate the coevolution of plants and herbivorous insects.

Its primary research objective is to establish a mechanistic framework connecting:

**Genotype → Functional Traits → Individual Interactions → Multi-species Interaction Networks**

The model aims to investigate:

1. The evolution of plant chemical defenses and volatile-mediated signals.
2. The evolution of insect sensory recognition and detoxification-related functional traits.
3. Reciprocal selection between plants and herbivorous insects.
4. The evolution of host preference, host use, and host specialization.
5. Host shifts and ecological divergence.
6. The effects of spatial structure and gene flow on evolutionary trajectories.
7. The emergence and reorganization of multispecies interaction networks.

## 1.2 Modeling Philosophy

EvoSim is a **theoretical, mechanism-driven Agent-Based Model (ABM)**.

It is not intended to be an empirical model fitted to one particular plant or insect species.

The model uses abstract genotypes and theoretical genotype-to-phenotype mappings to investigate general evolutionary mechanisms.

Therefore:

- Do not assume that every model trait must correspond to a directly measured biological gene.
- Do not force abstract model parameters to represent real physiological measurements.
- Do not introduce unnecessary molecular complexity.
- Do not add biological mechanisms merely to make the model appear more realistic.
- Every biological mechanism must have a clearly stated theoretical or empirical justification.
- Every simplification must have a clearly defined interpretation and scope.

The objective is scientific explanatory power, not maximum model complexity.

---

# 2. Scientific Integrity: Highest Priority

EvoSim is a scientific research project. Scientific correctness takes priority over code brevity, convenience, and implementation speed.

## 2.1 No Unauthorized Scientific Changes

Unless explicitly requested, Codex MUST NOT change:

- Biological assumptions.
- Mathematical equations.
- Functional relationships between traits.
- Genotype-to-phenotype mappings.
- Selection mechanisms.
- Fitness functions.
- Reproduction and inheritance rules.
- Mutation and recombination mechanisms.
- Population transition rules.
- Spatial dispersal mechanisms.
- Default parameter values.
- Parameter meanings or units.
- Life-cycle ordering.
- Scientific interpretation of output variables.

A proposed scientific improvement must first be documented separately with:

1. The current mechanism.
2. The biological or theoretical problem.
3. Supporting literature or mathematical reasoning.
4. The proposed alternative.
5. Expected consequences for model behavior.

Do not silently implement scientific changes during code refactoring.

## 2.2 Distinguish Code Refactoring from Model Development

Every task must be classified as one of the following:

- **Refactoring:** Improves code organization without changing scientific behavior.
- **Bug Fix:** Corrects an implementation error while preserving the intended scientific model.
- **Scientific Modification:** Changes a biological assumption, mathematical mechanism, or parameterization.
- **Analysis Extension:** Adds new methods for analyzing existing simulation outputs.

Do not combine these categories without explicit approval.

## 2.3 Preserve Existing Research Results

The existing v0.1.11 implementation and its associated outputs must be treated as an important historical baseline.

Never overwrite or delete:

- Historical simulation results.
- Historical configurations.
- Original input datasets.
- Previous model implementations.
- Archived scientific analyses.

When retiring source files or configurations, preserve a complete historical copy with:

- Original file paths.
- SHA256 checksums.
- Version information.
- Execution instructions.
- Necessary dependencies.

Historical results must remain independently interpretable.

---

# 3. Target Project Architecture

The following is the required long-term target structure.

New code must follow this architecture.

Existing files must be migrated gradually. Do not move all files simultaneously or create empty placeholder modules merely to match the directory tree.

## 3.1 Directory Structure

```text
EvoSim/
│
├── AGENTS.md
├── README.md
├── pyproject.toml
├── requirements.txt
├── example.py
│
├── docs/
│   ├── architecture.md
│   ├── model.md
│   ├── parameters.md
│   ├── results.md
│   └── development.md
│
├── evo_sim/
│   │
│   ├── core/
│   │   ├── __init__.py
│   │   ├── config.py
│   │   ├── controller.py
│   │   ├── initializer.py
│   │   ├── io.py
│   │   ├── registry.py
│   │   ├── runner.py
│   │   └── substep.py
│   │
│   └── models/
│       │
│       └── plant_insect/
│           │
│           ├── __init__.py
│           ├── simulator.py
│           ├── batch.py
│           │
│           ├── config/
│           │   ├── __init__.py
│           │   ├── loader.py
│           │   ├── schema.py
│           │   └── input_loader.py
│           │
│           ├── configs/
│           │   └── *.yaml
│           │
│           ├── inputs/
│           │   └── current/
│           │       └── *.csv
│           │
│           ├── biology/
│           │   ├── __init__.py
│           │   ├── agents.py
│           │   ├── traits.py
│           │   ├── spatial.py
│           │   │
│           │   ├── genetics/
│           │   │   ├── __init__.py
│           │   │   ├── qtl.py
│           │   │   ├── inheritance.py
│           │   │   ├── mutation.py
│           │   │   ├── recombination.py
│           │   │   ├── slim_backend.py
│           │   │   └── slim_reproduction.slim
│           │   │
│           │   ├── plants/
│           │   │   ├── __init__.py
│           │   │   ├── chemistry.py
│           │   │   ├── defense.py
│           │   │   └── demography.py
│           │   │
│           │   └── insects/
│           │       ├── __init__.py
│           │       ├── sensory.py
│           │       ├── host_selection.py
│           │       ├── feeding.py
│           │       └── demography.py
│           │
│           ├── simulation/
│           │   ├── __init__.py
│           │   ├── runtime.py
│           │   ├── random.py
│           │   ├── lifecycle.py
│           │   │
│           │   └── stages/
│           │       ├── __init__.py
│           │       ├── adult_reproduction.py
│           │       ├── oviposition.py
│           │       ├── offspring_genetics.py
│           │       ├── larval_development.py
│           │       ├── plant_dynamics.py
│           │       ├── plant_reproduction.py
│           │       └── recruitment.py
│           │
│           ├── records/
│           │   ├── __init__.py
│           │   ├── observations.py
│           │   ├── events.py
│           │   └── contributions.py
│           │
│           ├── audit/
│           │   ├── __init__.py
│           │   ├── run.py
│           │   ├── lifecycle.py
│           │   └── lifecycle_checks/
│           │       ├── __init__.py
│           │       ├── population.py
│           │       ├── genetics.py
│           │       ├── feeding.py
│           │       ├── resources.py
│           │       └── events.py
│           │
│           └── analysis/
│               ├── __init__.py
│               ├── common.py
│               ├── lineage_observations.py
│               │
│               ├── summaries/
│               │   ├── __init__.py
│               │   ├── lifecycle.py
│               │   ├── lineages.py
│               │   └── qtl.py
│               │
│               └── plots/
│                   ├── __init__.py
│                   └── diagnostics.py
│
├── tests/
│   ├── __init__.py
│   ├── __main__.py
│   ├── support.py
│   │
│   ├── fixtures/
│   │
│   ├── unit/
│   │   ├── test_config.py
│   │   ├── test_genetics.py
│   │   ├── test_chemistry.py
│   │   ├── test_host_selection.py
│   │   ├── test_feeding.py
│   │   ├── test_demography.py
│   │   └── test_spatial.py
│   │
│   ├── integration/
│   │   ├── test_runtime.py
│   │   ├── test_lifecycle.py
│   │   └── test_records.py
│   │
│   ├── scientific/
│   │   ├── test_scientific_boundaries.py
│   │   ├── test_mating.py
│   │   ├── test_plant_resources.py
│   │   └── test_independent_audit.py
│   │
│   └── acceptance/
│       ├── test_acceptance.py
│       └── test_batch.py
│
└── codex/
    ├── validation/
    └── history/
```

This structure is the architectural target, not permission to move or rewrite the entire repository immediately.

---

# 4. Module Responsibilities

Every directory has a clearly defined responsibility.

## 4.1 core/

Contains generic simulation infrastructure.

Responsibilities:

- Simulation scheduling.
- Generic controllers.
- Generic configuration interfaces.
- Generic input/output infrastructure.
- Generic substep interfaces.
- Model registration and execution.

Restrictions:

- Must not contain plant-specific or insect-specific biological equations.
- Must not import `plant_insect` modules.
- Must remain reusable by future EvoSim models.

## 4.2 config/

Responsible for:

- Loading YAML configurations.
- Validating configuration schemas.
- Loading input CSV files.
- Checking configuration compatibility.
- Resolving configuration inheritance.

Do not implement biological equations here.

Keep experiment configurations in `configs/` and biological input data in `inputs/`.

## 4.3 biology/

Contains the actual scientific mechanisms of the model.

This directory answers:

**What biological processes does EvoSim represent, and how are they mathematically described?**

### biology/genetics/

Responsible for:

- QTL architecture.
- Genotype representation.
- Allele inheritance.
- Recombination.
- Mutation.
- SLiM integration.
- Genotype transmission.

### biology/plants/

Responsible for:

- Plant chemical traits.
- Chemical defense.
- Resource-related plant traits.
- Plant survival and demographic mechanisms.

### biology/insects/

Responsible for:

- Insect sensory traits.
- Host recognition.
- Host selection and rejection.
- Feeding and resource utilization.
- Insect survival and demography.

### biology/traits.py

Responsible for the mapping between genetic information and expressed functional traits.

Do not duplicate genotype-to-phenotype calculations in different biological modules.

### biology/spatial.py

Responsible for:

- Spatial representation.
- Dispersal kernels.
- Spatial neighborhood calculations.
- Spatial constraints on interactions.

Do not implement life-cycle scheduling here.

## 4.4 simulation/

Responsible for executing the model.

The simulation layer organizes the sequence of biological events but does not redefine their scientific equations.

### simulation/lifecycle.py

The central generation controller.

It calls individual life-cycle stages in the correct order.

### simulation/stages/

Each file represents one clearly defined stage.

The existing scientific life-cycle order must be preserved during refactoring:

1. Adult mating and reproduction.
2. Per-egg host selection and oviposition.
3. Offspring genetic inheritance through SLiM.
4. Larval feeding and survival.
5. Plant survival and demographic updates.
6. Plant reproduction and resource recovery.
7. Recruitment of the next generation.

The actual existing implementation must be inspected before assigning code to these stages.

Do not change the order simply because a different order appears cleaner.

## 4.5 records/

Responsible for recording simulation outputs.

Includes:

- Population observations.
- Individual-level events.
- Reproductive events.
- Genetic contributions.
- Stage-level records.

Recording code must not independently calculate biological outcomes.

## 4.6 audit/

Responsible for independently checking model execution.

This directory is scientifically important.

The audit system must independently reconstruct or verify critical model behavior using recorded data and configuration information.

**Do not make the audit call the production implementation's biological calculation functions to verify those same functions.**

Otherwise, the same error could occur in both the simulator and its audit, causing a false successful validation.

Preserve the independence of `audit_lifecycle.py` when splitting it into smaller modules.

## 4.7 analysis/

Responsible for interpreting simulation outputs.

Includes:

- Data summaries.
- Lineage observations.
- Genetic summaries.
- Interaction statistics.
- Visualization.

Restrictions:

- Must not modify simulation state.
- Must not implement biological mechanisms.
- Must not silently change raw simulation records.
- Must preserve historical output compatibility.

## 4.8 tests/

Tests must be organized by their purpose.

- `unit/`: individual functions and mathematical mechanisms.
- `integration/`: interactions between modules.
- `scientific/`: biological and mathematical invariants.
- `acceptance/`: complete simulation and batch behavior.
- `fixtures/`: reusable test inputs and reference data.

Keep the existing test execution interface working.

The command `python -m tests` must continue to run the complete test suite unless explicitly approved otherwise.

---

# 5. Module Dependency Rules

The dependency direction must remain clear.

The following rules are mandatory:

1. `core/` must not depend on biological model modules.
2. `biology/` must not depend on `simulation/`, `analysis/`, or `audit/`.
3. `simulation/` may call `biology/` functions to execute model processes.
4. `records/` may receive outputs from `simulation/`, but must not redefine biological mechanisms.
5. `audit/` may read configurations and records, but must independently verify critical calculations.
6. `analysis/` may read simulation outputs and records, but must not modify the original data.
7. `config/` must not depend on analysis or plotting modules.
8. Avoid circular imports.
9. Avoid hidden dependencies through mutable global variables.
10. Pass important state and dependencies explicitly.

If a dependency violates these rules, explain why before introducing it.

---

# 6. Code Style and Naming Conventions

## 6.1 Language

- Python code: English.
- Function and variable names: English.
- Class names: English.
- Scientific equations and terminology: use standard English terminology.
- Documentation: Chinese or English, depending on the intended audience.
- Comments explaining biological assumptions may be written in Chinese for researcher accessibility.

Use UTF-8.

## 6.2 Naming

Follow standard Python conventions:

- Variables and functions: `snake_case`.
- Classes: `PascalCase`.
- Constants: `UPPER_SNAKE_CASE`.
- Modules and packages: lowercase.

Use descriptive names.

Avoid unexplained abbreviations.

Prefer:

`calculate_larval_survival_probability()`

Instead of:

`calc_lsurv()`

Prefer:

`plant_resource_deficit`

Instead of:

`prd`

Biological abbreviations such as QTL, VOC, PSM, OR, OBP, P450, GST, and UGT are acceptable.

## 6.3 Function Design

Each function should have one clearly defined responsibility.

Avoid:

- Very long functions.
- Excessively nested conditionals.
- Functions that combine unrelated biological processes.
- Hidden changes to global state.
- Large numbers of positional arguments.
- Unnecessary one-line helper functions.
- Excessive abstraction.

Do not fragment code merely to reduce line count.

The goal is conceptual clarity.

## 6.4 Type Hints

Use type annotations for important functions and public interfaces.

Example:

```python
def calculate_host_utility(
    insect_traits: InsectTraits,
    plant_traits: PlantTraits,
) -> float:
    """Calculate the utility of a candidate host plant."""
```

Use dataclasses or typed structures where they genuinely improve state representation.

Do not introduce complex class hierarchies without a clear need.

## 6.5 Docstrings

Every important function must have a concise docstring.

A scientific function should describe:

1. What biological process it represents.
2. What mathematical operation it performs.
3. Input variables and units.
4. Output variables and units.
5. Relevant assumptions or limitations.

Example:

```python
def calculate_defense_processing_match(
    insect_processing: np.ndarray,
    plant_defense: np.ndarray,
) -> float:
    """
    Calculate the functional match between insect processing capacity
    and plant chemical defense.

    The current model uses a limiting-channel formulation:

        q = min(1, min_k(I_k / D_k))

    Parameters
    ----------
    insect_processing:
        Abstract insect processing capacities.
    plant_defense:
        Abstract plant defense burdens.

    Returns
    -------
    float
        Relative feeding utilization capacity.

    Notes
    -----
    This is an abstract functional matching mechanism, not a direct
    measurement of enzyme activity or detoxification rate.
    """
```

## 6.6 Comments

Comments should explain:

- Why a biological assumption is used.
- Why a mathematical relationship is selected.
- Why a particular implementation decision is necessary.
- What limitations exist.

Do not write comments that merely repeat the code.

---

# 7. Scientific Equations and Parameter Management

## 7.1 Keep Equations Identifiable

Important scientific equations must be easy to locate.

Do not hide important equations inside long nested expressions.

Use intermediate variables with meaningful names.

For example:

```python
resource_deficit = 1.0 - post_feeding_resource / resource_capacity

mortality_hazard = (
    background_mortality_rate
    + resource_deficit_mortality_rate * resource_deficit
)

survival_probability = np.exp(-mortality_hazard)
```

The implementation must preserve the existing equation unless a scientific modification has been approved.

## 7.2 Parameter Rules

Never introduce arbitrary numerical values directly into biological functions.

Use configuration parameters or explicitly documented mathematical constants.

Every scientific parameter should have:

- Name.
- Biological or mathematical meaning.
- Unit or dimensionless status.
- Current value.
- Source or theoretical justification.
- Status of evidence.
- Known limitations.

Distinguish:

- Empirically calibrated parameters.
- Literature-informed parameters.
- Theoretical parameters.
- Implementation parameters.
- Observation and analysis settings.

A literature citation does not automatically validate the numerical value used in EvoSim.

Do not invent biological evidence.

If evidence is insufficient, label the parameter as a theoretical assumption.

## 7.3 Parameter Changes

Changing a default parameter requires:

1. A scientific justification.
2. An update to `docs/parameters.md`.
3. An update to the relevant configuration.
4. Appropriate tests.
5. A clear description of the expected model consequences.

Never silently change default parameter values during refactoring.

---

# 8. Randomness and Reproducibility

Randomness is part of the scientific model.

The following must be preserved during behavior-preserving refactoring:

- Random seed derivation.
- Random number generator type.
- Random number generator streams.
- Order of random draws.
- Random sampling algorithms.
- SLiM random behavior.
- Reproducibility of fixed-seed runs.

Do not replace one random sampling function with another merely because the replacement looks cleaner.

Do not move random draws across life-cycle stages.

Do not introduce additional random draws during refactoring.

If a random-number implementation must change, treat it as a separate scientific and computational change.

---

# 9. Testing and Validation

Testing must cover both software correctness and scientific correctness.

## 9.1 Before Refactoring

Before changing an existing module:

1. Inspect its implementation.
2. Inspect all relevant tests.
3. Identify its scientific equations.
4. Identify its inputs and outputs.
5. Record baseline behavior.
6. Use fixed random seeds where applicable.

## 9.2 During Refactoring

Refactor one logical module at a time.

After each change:

1. Run relevant unit tests.
2. Compare outputs with the original implementation.
3. Compare intermediate records where necessary.
4. Check random-number behavior.
5. Check SLiM integration if affected.
6. Review all changed imports and interfaces.

## 9.3 Final Validation

Run:

```bash
python -m ruff check evo_sim tests example.py

python -m ruff format --check evo_sim tests example.py

python -m tests
```

Install development dependencies if necessary:

```bash
python -m pip install -e ".[dev]"
```

Do not claim a test passed unless it was actually executed.

Distinguish:

- Not run.
- Passed.
- Failed.
- Blocked by missing dependencies.

## 9.4 Scientific Audit

A successful unit-test suite does not prove that the scientific model is biologically valid.

Keep separate:

- Software implementation validation.
- Mathematical consistency checks.
- Scientific mechanism validation.
- Parameter uncertainty analysis.
- Empirical validation.

Do not claim scientific validation based only on successful software tests.

---

# 10. Refactoring Workflow

The repository must be refactored incrementally.

Do not perform a full-project rewrite in one operation.

## Phase 1: Architecture Inspection

No code changes during this phase.

Tasks:

1. Read the existing repository.
2. Map all current modules.
3. Identify module responsibilities.
4. Identify dependencies.
5. Identify long and overloaded functions.
6. Identify duplicate calculations.
7. Compare the existing structure with the target architecture.
8. Produce `docs/architecture.md`.
9. Produce a detailed migration plan.

The report must distinguish actual problems from optional improvements.

## Phase 2: Establish the New Architecture

Create the necessary package structure.

Update:

- `AGENTS.md`.
- `docs/architecture.md`.
- Existing import and package conventions.

Do not move large groups of files before the migration plan is approved.

## Phase 3: Refactor Biological Mechanisms

Recommended order:

1. Genetics.
2. Plant chemistry and defense.
3. Insect sensory and host selection.
4. Feeding and resource utilization.
5. Plant and insect demography.
6. Spatial mechanisms.

Each module must have clear biological responsibilities.

## Phase 4: Refactor Simulation Lifecycle

Refactor the existing lifecycle controller into independent stages.

Preserve the exact existing execution order.

Keep the central lifecycle controller concise and readable.

## Phase 5: Refactor Records, Audit, and Analysis

Preserve:

- Original output schemas.
- Historical filenames.
- Event identifiers.
- Individual identifiers.
- Independent audit logic.
- Existing scientific interpretation.

## Phase 6: Tests and Documentation

Update test locations and imports.

Update all relevant documentation and execution instructions.

Run the complete validation suite.

Do not release a new scientific model version merely because the source code has been reorganized.

---

# 11. File Migration Rules

The following mapping provides the initial migration guide.

| Existing file | Target location |
|---|---|
| `helpers/agents.py` | `biology/agents.py` |
| `helpers/genetics.py` | `biology/genetics/` |
| `helpers/slim_backend.py` | `biology/genetics/slim_backend.py` |
| `helpers/slim_reproduction.slim` | `biology/genetics/slim_reproduction.slim` |
| `helpers/chemistry.py` | `biology/plants/chemistry.py` |
| `helpers/voc.py` | Determine the correct plant signal or insect sensory module after inspecting its actual responsibility |
| `helpers/ecology.py` | Split into relevant host selection and feeding modules |
| `helpers/demography.py` | Split into plant and insect demography modules |
| `helpers/spatial.py` | `biology/spatial.py` |
| `helpers/traits.py` | `biology/traits.py` |
| `helpers/modules.py` | Appropriate functional trait module after reviewing its responsibilities |
| `helpers/config.py` | `config/loader.py` and `config/schema.py` |
| `helpers/inputs.py` | `config/input_loader.py` |
| `helpers/runtime.py` | `simulation/runtime.py` |
| `helpers/observations.py` | `records/observations.py` |
| `helpers/contributions.py` | `records/contributions.py` |
| `substeps/lifecycle.py` | `simulation/lifecycle.py` and `simulation/stages/` |
| `audit.py` | `audit/run.py` |
| `audit_lifecycle.py` | `audit/lifecycle.py` and independent check modules |
| `analysis/summarize_lifecycle.py` | `analysis/summaries/lifecycle.py` |
| `analysis/summarize_lineages.py` | `analysis/summaries/lineages.py` |
| `analysis/summarize_qtl_run.py` | `analysis/summaries/qtl.py` |
| `analysis/plot_diagnostics.py` | `analysis/plots/diagnostics.py` |

This is a responsibility mapping, not an instruction to blindly move files.

Inspect the actual functions before splitting them.

Preserve public interfaces or provide compatibility wrappers when necessary.

Do not create unnecessary modules.

Do not create new generic folders named `helpers/`, `utils/`, or `misc/`.

---

# 12. Data, Results, and Historical Archives

## 12.1 Research Results

`results/` is reserved for formal scientific experiments explicitly requested by the user.

Development and validation runs must not be written into `results/`.

Use:

```text
codex/validation/v<version>/run<N>/
```

for development runs, smoke tests, short mechanism checks, and implementation validation.

Directories must not overwrite previous runs.

## 12.2 Historical Archives

Use:

```text
codex/history/
```

for preserved historical source files and development records.

The `codex/` directory is Git-ignored.

Important scientific documentation must remain in `docs/`, not only in local Codex records.

Do not automatically delete or migrate existing research results.

---

# 13. Git and Version Control

After completing and validating a coherent batch of changes:

1. Review the complete Git diff.
2. Ensure only relevant files are included.
3. Commit the changes.
4. Write a descriptive commit message.
5. Record the actual tests executed.
6. State any remaining limitations.

Do not create a formal release tag for an unvalidated scientific model.

Formal versions must use:

```text
v<version>
```

Annotated tags must point to the corresponding committed source code.

Package version, configuration protocol, and documentation must remain consistent.

Never rewrite shared commits or tags merely to update descriptions.

At the end of a development task, report:

- Files changed.
- Main modifications.
- Scientific behavior preserved or changed.
- Tests actually executed.
- Commit hash.
- Tag status, if applicable.
- Remote synchronization status.

---

# 14. Documentation Responsibilities

Keep existing scientific documentation up to date.

| File | Responsibility |
|---|---|
| `README.md` | Project overview, installation, execution |
| `docs/architecture.md` | Package architecture, module responsibilities, dependencies |
| `docs/model.md` | Scientific model, equations, biological assumptions |
| `docs/parameters.md` | Parameter definitions, values, units, evidence |
| `docs/results.md` | Output formats, analysis, plotting, commands |
| `docs/development.md` | Development workflow, testing, refactoring |

Do not create duplicate versioned documentation without a specific reason.

When a scientific mechanism changes, update `docs/model.md`.

When a parameter changes, update `docs/parameters.md`.

When the code architecture changes, update `docs/architecture.md`.

When output schemas or analysis commands change, update `docs/results.md`.

---

# 15. Codex Working Protocol

For every non-trivial development task, follow this process:

## Step 1: Understand

Read the relevant source code, configuration, documentation, and tests.

Do not assume that filenames accurately describe their current responsibilities.

## Step 2: Explain

Before making substantial changes, briefly explain:

- What is wrong.
- Why it matters.
- Which files are involved.
- What will change.
- What will explicitly remain unchanged.

## Step 3: Plan

Divide large tasks into small, independently reviewable changes.

## Step 4: Implement

Modify only the files necessary for the current task.

Avoid unrelated cleanup.

## Step 5: Validate

Run appropriate tests and inspect the resulting differences.

## Step 6: Report

Clearly state what was changed, what was tested, and what remains incomplete.

Never report an unexecuted task as completed.

---

# 16. Final Development Principles

All development must follow these principles:

1. Scientific correctness over coding convenience.
2. Clear biological meaning over excessive abstraction.
3. Readability over minimum line count.
4. Explicit dependencies over hidden state.
5. Reproducibility over implementation shortcuts.
6. Independent auditing over self-confirming tests.
7. Incremental refactoring over full-project rewrites.
8. Documented assumptions over unexplained numerical values.
9. Stable interfaces over unnecessary API changes.
10. Long-term scientific maintainability over short-term coding speed.

**EvoSim must remain understandable to a biological researcher who wants to understand the evolutionary mechanisms, not only to a software engineer who wants to execute the code.**

The directory structure, function names, documentation, and tests should collectively make the scientific model transparent.

---

# 17. Architecture Migration Status (2026-11-06)

The approved behavior-preserving migration uses the packages documented in
`docs/architecture.md`. Biological implementations now live in `biology/`,
ordered stages in `simulation/stages/`, and independent checks in
`audit/lifecycle_checks/`. Tests are grouped by purpose and recursively discovered
by `python -m tests`. The v0.1.11 scientific protocol and original release tag
remain the baseline. Keep the thin historical analysis command entry points;
they forward to current implementations and do not duplicate scientific logic.


# 18. Scientific Protocol v0.1.12 (2026-11-09)

The user explicitly authorized the ecological scientific modification documented in
`docs/ecology_v0.1.12_design.md`. Current plant coordinates are PSM1, PSM2, Z,
VOC1, VOC2; insect coordinates have ten axes including two raw contact responses.
Processing acts after ingestion, never on food demand. Defoliation is actual
consumption / pre-feeding stock; tolerance changes recovery, and recovered stock
supports reproduction. Retain independent audit reconstruction and transactional
publication. Preserve v0.1.11 commit fd9ba73, its inputs and historical semantics;
new readers must not reinterpret its outputs. SLiM's core genetic algorithm and
speciation/reproductive-isolation mechanisms remain unchanged.
