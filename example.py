"""Run the current plant–insect model through the minimal Runner API."""

from evo_sim.models.plant_insect import get_runner

if __name__ == "__main__":
    runner = get_runner()
    runner.init()
    runner.step()
    print(runner.result["metadata"]["status"])
