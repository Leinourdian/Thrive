﻿using System.Collections.Generic;
using AutoEvo;

/// <summary>
///   Custom auto-evo run for use in editor simulations to predict population numbers given specified changes
/// </summary>
public class EditorAutoEvoRun : AutoEvoRun
{
    public EditorAutoEvoRun(GameWorld world, Species originalEditedSpecies, Species modifiedProperties) : base(world)
    {
        OriginalEditedSpecies = originalEditedSpecies;
        ModifiedProperties = modifiedProperties;
    }

    public Species OriginalEditedSpecies { get; }
    public Species ModifiedProperties { get; }

    protected override void GatherInfo(Queue<IRunStep> steps)
    {
        // Custom run setup for editor's use
        var map = Parameters.World.Map;

        ModifiedProperties.Mutated = true;

        steps.Enqueue(new CalculatePopulation(SimulationParameters.Instance.AutoEvoConfiguration, map,
            new List<Species> { ModifiedProperties },
            new List<Species> { OriginalEditedSpecies }, true) { CanRunConcurrently = false });

        //if (ModifiedProperties.Population <= 0)
        //{
        //    ModifiedProperties.Population = 2;
        //}

        AddPlayerSpeciesPopulationChangeClampStep(steps, map,
            OriginalEditedSpecies.PlayerSpecies ? ModifiedProperties : null, OriginalEditedSpecies);
    }
}
