namespace AutoEvo
{
    using System;

    public abstract class FoodSource
    {
        private readonly Compound glucose = SimulationParameters.Instance.GetCompound("glucose");
        private readonly Compound atp = SimulationParameters.Instance.GetCompound("atp");

        public abstract float TotalEnergyAvailable();

        /// <summary>
        ///   Provides a fitness metric to determine population adjustments for species in a patch.
        /// </summary>
        /// <param name="microbe">The species to be evaluated.</param>
        /// <param name="simulationCache">
        ///   Cache that should be used to reduce amount of times expensive computations are run
        /// </param>
        /// <returns>
        ///   A float to represent score. Scores are only compared against other scores from the same FoodSource,
        ///   so different implementations do not need to worry about scale.
        /// </returns>
        public abstract (float, float) FitnessScore(Species microbe, SimulationCache simulationCache);

        /// <summary>
        ///   A description of this niche. Needs to support translations changing and be player readable
        /// </summary>
        /// <returns>A formattable that has the description in it</returns>
        public abstract IFormattable GetDescription();

        // changed: added cache, patch, bool
        protected (float, float) EnergyGenerationScore(MicrobeSpecies species, Compound compound, SimulationCache cache, Patch patch,
            bool compoundIsSearchable, float compoundAmount)
        {
            // add floor calculation
            float glucoseCreated = 0.0f;
            float glucoseInput = 0.0f;
            float glucoseATPOutput = 0.0f; // silly name, it's actually atpOutputFromGlucose
            float compoundInput = 0.0f;
            float compoundATPOutput = 0.0f;
            var energyCreationScore = 0.0f;

            foreach (var organelle in species.Organelles)
            {
                foreach (var process in organelle.Definition.RunnableProcesses)
                {
                    if (process.Process.Inputs.ContainsKey(glucose) && compound != glucose)
                    {
                        glucoseInput += process.Process.Inputs[glucose];
                        glucoseATPOutput += process.Process.Outputs[atp];
                    }
                    else if (process.Process.Inputs.ContainsKey(compound))
                    {
                        compoundInput += process.Process.Inputs[compound];

                        if (process.Process.Outputs.ContainsKey(glucose))
                        {
                            glucoseCreated += process.Process.Outputs[glucose];
                        }

                        if (process.Process.Outputs.ContainsKey(atp)) //TryGetValue(atp, out var atpAmount))
                        {
                            compoundATPOutput += process.Process.Outputs[atp];
                        }
                    }
                }
            }

            //float glucoseAmount = 0.0f;
            //if (compoundInput != 0.0f)
            //{
            //    glucoseAmount += compoundAmount * glucoseCreated / compoundInput;
            //}

            //if (glucoseInput != 0.0f)
            //{
            //    energyCreationScore += glucoseAmount * glucoseATPOutput / glucoseInput;
            //}

            //if (compoundATPOutput != 0.0f)
            //{
            //    energyCreationScore += compoundAmount * compoundATPOutput / compoundInput;
            //}

            float glucoseAmount = 0.0f;
            if (compoundInput != 0.0f)
            {
                glucoseAmount += compoundAmount * glucoseCreated / compoundInput;
            }

            if (glucoseInput != 0.0f)
            {
                energyCreationScore += Math.Min(glucoseCreated, glucoseInput) * glucoseATPOutput / glucoseInput;
            }

            if (compoundATPOutput != 0.0f)
            {
                energyCreationScore += compoundATPOutput;
            }

            //if (energyCreationScore >= 0.0f)
            //{
            //}
            //else
            //{
            //    Console.WriteLine("glucoseAmount: " + glucoseAmount +
            //        "\nglucoseCreated: " + glucoseCreated +
            //        "\ncompoundAmount: " + compoundAmount +
            //        "\ncompoundInput: " + compoundInput +
            //        "\ncompoundOutput: " + compoundOutput +
            //        "\nglucoseInput: " + glucoseInput +
            //        "\nglucoseOutput: " + glucoseOutput +
            //        "\nenergyCreationScore: " + energyCreationScore +
            //        "\nspecies: " + species.FormattedName);
            //}

            return (energyCreationScore, compoundInput);
            //return (energyCreationScore, compoundATPOutput, compoundATPOutput / compoundInput);



            //var energyCreationScore = 0.0f;
            //foreach (var organelle in species.Organelles)
            //{
            //    foreach (var process in organelle.Definition.RunnableProcesses)
            //    {
            //        if (process.Process.Inputs.ContainsKey(compound))
            //        {
            //            if (process.Process.Outputs.ContainsKey(glucose))
            //            {
            //                energyCreationScore += process.Process.Outputs[glucose]
            //                    / process.Process.Inputs[compound] * Constants.AUTO_EVO_GLUCOSE_USE_SCORE_MULTIPLIER;
            //            }

            //            if (process.Process.Outputs.ContainsKey(atp))
            //            {
            //                energyCreationScore += process.Process.Outputs[atp]
            //                    / process.Process.Inputs[compound] * Constants.AUTO_EVO_ATP_USE_SCORE_MULTIPLIER;
            //            }
            //        }
            //    }
            //}

            //// changed: added following to make this about collecting.
            ////          not yet but balance this better.
            ////          should the minimum be higher than 0?
            //if (energyCreationScore > 0 && compoundIsSearchable)
            //{
            //    var energyBalanceInfo = cache.GetEnergyBalanceForSpecies(species, patch);
            //    float actualSpeed = species.BaseSpeed * Math.Min(1.0f, energyBalanceInfo.FinalBalanceStationary /
            //        (energyBalanceInfo.TotalConsumption - energyBalanceInfo.TotalConsumptionStationary));
            //    float result = Math.Max(0, species.BaseHexSize * species.BaseHexSize * actualSpeed * actualSpeed);
            //    return result;
            //}

            //return energyCreationScore;
        }
    }
}
