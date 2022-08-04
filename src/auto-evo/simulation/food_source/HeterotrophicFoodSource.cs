﻿namespace AutoEvo
{
    using System;

    public class HeterotrophicFoodSource : FoodSource
    {
        private readonly Compound oxytoxy = SimulationParameters.Instance.GetCompound("oxytoxy");

        private readonly MicrobeSpecies prey;
        private readonly Patch patch;
        private readonly float preyHexSize;
        private readonly float preySpeed;
        private readonly float totalEnergy;

        public HeterotrophicFoodSource(Patch patch, MicrobeSpecies prey, SimulationCache simulationCache)
        {
            this.prey = prey;
            this.patch = patch;
            preyHexSize = simulationCache.GetBaseHexSizeForSpecies(prey);
            preySpeed = simulationCache.GetBaseSpeedForSpecies(prey);
            patch.SpeciesInPatch.TryGetValue(prey, out long population);
            // changed: not yet but this scales to billions.
            //          size needs to concider population calculation.
            //          this should be dependent on how much the species is preyed upon.
            //          preying should reduce prey population.
            //          this should use energyBySpecies to avoid exponential growth
            totalEnergy = population * prey.Organelles.Count * Constants.AUTO_EVO_PREDATION_ENERGY_MULTIPLIER;
        }

        public override (float, float) FitnessScore(Species species, SimulationCache simulationCache)
        {
            var microbeSpecies = (MicrobeSpecies)species;

            // No cannibalism
            if (microbeSpecies == prey)
            {
                return (0.0f, 0.0f);
            }

            var behaviourScore = microbeSpecies.Behaviour.Aggression / Constants.MAX_SPECIES_AGGRESSION;

            var microbeSpeciesHexSize = microbeSpecies.BaseHexSize;
            var predatorSpeed = microbeSpecies.BaseSpeed;
            predatorSpeed += simulationCache.GetEnergyBalanceForSpecies(microbeSpecies, patch.Biome).FinalBalance; // changed: not yet but ad to this

            // It's great if you can engulf this prey, but only if you can catch it
            var engulfScore = 0.0f;
            if (microbeSpeciesHexSize / preyHexSize >
                Constants.ENGULF_SIZE_RATIO_REQ && !microbeSpecies.MembraneType.CellWall) // changed: not yet but change to CanTryToEat
            {
                engulfScore = Constants.AUTO_EVO_ENGULF_PREDATION_SCORE;
            }

            engulfScore *= predatorSpeed > preySpeed ? 1.0f : Constants.AUTO_EVO_ENGULF_LUCKY_CATCH_PROBABILITY; // changed: not yet but behavior could matter

            var pilusScore = 0.0f; // changed: not yet but store this in species memory or smth to save time
            var oxytoxyScore = 0.0f;
            foreach (var organelle in microbeSpecies.Organelles)
            {
                if (organelle.Definition.HasPilusComponent)
                {
                    pilusScore += Constants.AUTO_EVO_PILUS_PREDATION_SCORE;
                    continue;
                }

                foreach (var process in organelle.Definition.RunnableProcesses)
                {
                    if (process.Process.Outputs.TryGetValue(oxytoxy, out var oxytoxyAmount))
                    {
                        oxytoxyScore += oxytoxyAmount * Constants.AUTO_EVO_TOXIN_PREDATION_SCORE;
                    }
                }
            }

            // Pili are much more useful if the microbe can close to melee
            pilusScore *= predatorSpeed > preySpeed ? 1.0f : Constants.AUTO_EVO_ENGULF_LUCKY_CATCH_PROBABILITY;

            // predators are less likely to use toxin against larger prey, unless they are opportunistic
            if (preyHexSize > microbeSpeciesHexSize)
            {
                oxytoxyScore *= microbeSpecies.Behaviour.Opportunism / Constants.MAX_SPECIES_OPPORTUNISM;  // changed: not yet but could remove this
            }

            // Intentionally don't penalize for osmoregulation cost to encourage larger monsters
            return (0.0f, behaviourScore * (pilusScore + engulfScore + microbeSpeciesHexSize + oxytoxyScore)); // changed: not yet but debatable; it sounds fun though
        }

        public override IFormattable GetDescription()
        {
            return new LocalizedString("PREDATION_FOOD_SOURCE", prey.FormattedNameBbCode);
        }

        public override float TotalEnergyAvailable()
        {
            return totalEnergy;
        }
    }
}
