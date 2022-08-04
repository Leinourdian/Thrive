namespace AutoEvo
{
    using System;
    using Godot;

    public class CompoundFoodSource : FoodSource
    {
        private readonly Patch patch;
        private readonly Compound compound;
        private readonly float totalCompound;
        private readonly float compoundDensity;

        public CompoundFoodSource(Patch patch, Compound compound)
        {
            this.patch = patch;
            this.compound = compound;
            if (patch.Biome.Compounds.TryGetValue(compound, out var compoundData))
            {
                compoundDensity = compoundData.Density;
                totalCompound = compoundDensity * compoundData.Amount;
            }
            else
            {
                totalCompound = 0.0f;
            }
        }

        public override (float, float) FitnessScore(Species species, SimulationCache simulationCache)
        {
            var microbeSpecies = (MicrobeSpecies)species;
            float energyCoefficient = 100.0f;
            float effectiveDensity = 10.0f * totalCompound; //is this needed?

            // changed: not yet but make this about collecting compounds
            //          collectionScore = speed * size ( and maybe fighting ability too).
            //          added parameters cache and patch
            //var compoundUseScore = EnergyGenerationScore(microbeSpecies, compound, simulationCache, patch, true);

            //var energyGenerated = EnergyGenerationScore(microbeSpecies, compound, simulationCache, patch, true, totalCompound);

            var energyBalanceInfo = simulationCache.GetEnergyBalanceForSpecies(microbeSpecies, patch.Biome);

            float actualSpeed = microbeSpecies.BaseSpeed * Mathf.Clamp(energyBalanceInfo.FinalBalanceStationary /
                (energyBalanceInfo.TotalConsumption - energyBalanceInfo.TotalConsumptionStationary), 0.0f, 1.0f);

            float collectingLevel = Mathf.Pow(microbeSpecies.BaseHexSize, 1.0f) * actualSpeed;
            float collectedCompound = collectingLevel * compoundDensity;

            // use this?
            float trueEnergyNeed = energyBalanceInfo.TotalConsumption + energyBalanceInfo.FinalBalance;

            (float energyGenerated, float compoundUsed) = EnergyGenerationScore(microbeSpecies, compound, simulationCache,
                patch, true, collectedCompound);

            if (energyGenerated <= 0.0f)
                return (0.0f, 1.0f);

            float collectedEnergy = collectedCompound * energyGenerated / compoundUsed;

            float energyRatio = trueEnergyNeed / energyGenerated;
            float compoundNeedPerSecond = energyRatio * compoundUsed;
            float compoundCollectedPerSecond = effectiveDensity / collectingLevel;

            // we wanna figure out how much stuff you get at a certain density
            // calculate energy need, floor is when collectedEnergy * currentDensity = energy need
            // Should this use TotalConsumption instead? Debatable... you can still sometimes hunt at an energy deficit
            //float floor = Math.Min(1.0f, energyBalanceInfo.TotalConsumptionStationary / 100);
            //float floor = Math.Min(1.0f, energyBalanceInfo.TotalConsumptionStationary / energyPerSecond);
            //float floor = Math.Min(1.0f, energyCoefficient * compoundNeedPerSecond / actualSpeed);//actualSpeed);
            float floor = Math.Min(1.0f, microbeSpecies.BaseHexSize / 50);

            //Math.Min(1.0f, energyBalanceInfo.TotalConsumptionStationary / collectingLevel);

            // May wanna compromise later
            if (trueEnergyNeed > energyGenerated || floor >= 1.0f)
                return (0.0f,  1.0f);

            //return (collectedEnergy, floor);
            return (collectingLevel, floor);







            //var energyCost = simulationCache
            //    .GetEnergyBalanceForSpecies(microbeSpecies, patch)
            //    .TotalConsumptionStationary; // changed: not yet but should this take movement cost into account?

            //return compoundUseScore;// / energyCost; // changed: not yet but is this good?
        }

        public override IFormattable GetDescription()
        {
            // TODO: somehow allow the compound name to translate properly. Maybe we need to use bbcode to refer to the
            // compounds?
            return new LocalizedString("COMPOUND_FOOD_SOURCE", compound.Name);
        }

        public override float TotalEnergyAvailable()
        {
            return totalCompound * Constants.AUTO_EVO_COMPOUND_ENERGY_AMOUNT;
        }
    }
}
