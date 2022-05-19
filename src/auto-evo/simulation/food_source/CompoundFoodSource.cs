namespace AutoEvo
{
    using System;

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
                totalCompound = 1.0f; // 0.0f;
            }
        }

        public override (float, float) FitnessScore(Species species, SimulationCache simulationCache)
        {
            var microbeSpecies = (MicrobeSpecies)species;
            float energyCoefficient = 1.0f; //this doesn't seem to matter at all

            // changed: not yet but make this about collecting compounds
            //          collectionScore = speed * size ( and maybe fighting ability too).
            //          added parameters cache and patch
            //var compoundUseScore = EnergyGenerationScore(microbeSpecies, compound, simulationCache, patch, true);

            var energyGenerated = EnergyGenerationScore(microbeSpecies, compound, simulationCache, patch, true, totalCompound);

            var energyBalanceInfo = simulationCache.GetEnergyBalanceForSpecies(microbeSpecies, patch);

            float actualSpeed = microbeSpecies.BaseSpeed * Math.Max(0, Math.Min(1.0f, energyBalanceInfo.FinalBalanceStationary /
                (energyBalanceInfo.TotalConsumption - energyBalanceInfo.TotalConsumptionStationary)));

            float collectingLevel = (float)Math.Sqrt(microbeSpecies.BaseHexSize) * actualSpeed;
            float collectedCompound = 1.0f; //collectingLevel * compoundDensity;

            // should this use the constant multiplier?
            float collectedEnergy = energyCoefficient * EnergyGenerationScore(microbeSpecies, compound, simulationCache, patch, true, collectedCompound);

            // we wanna figure out how much stuff you get at a certain density
            // calculate energy need, floor is when collectedEnergy * currentDensity = energy need
            float floor = Math.Min(1.0f, energyBalanceInfo.TotalConsumptionStationary / collectedEnergy);




            return (collectedEnergy, floor);







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
