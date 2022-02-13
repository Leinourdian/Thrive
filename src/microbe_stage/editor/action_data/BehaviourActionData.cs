﻿using System;

[JSONAlwaysDynamicType]
public class BehaviourActionData : MicrobeEditorCombinableActionData
{
    public float OldValue;
    public float NewValue;
    public BehaviouralValueType Type;

    public BehaviourActionData(float oldValue, float newValue, BehaviouralValueType type)
    {
        OldValue = oldValue;
        NewValue = newValue;
        Type = type;
    }

    public override MicrobeActionInterferenceMode GetInterferenceModeWith(CombinableActionData other)
    {
        if (other is BehaviourActionData behaviourChangeActionData && behaviourChangeActionData.Type == Type)
        {
            // If the value has been changed back to a previous value
            if (Math.Abs(NewValue - behaviourChangeActionData.OldValue) < MathUtils.EPSILON &&
                Math.Abs(behaviourChangeActionData.NewValue - OldValue) < MathUtils.EPSILON)
                return MicrobeActionInterferenceMode.CancelsOut;

            // If the value has been changed twice
            if (Math.Abs(NewValue - behaviourChangeActionData.OldValue) < MathUtils.EPSILON ||
                Math.Abs(behaviourChangeActionData.NewValue - OldValue) < MathUtils.EPSILON)
                return MicrobeActionInterferenceMode.Combinable;
        }

        return MicrobeActionInterferenceMode.NoInterference;
    }

    public override int CalculateCost()
    {
        // TODO: should this be free?
        return 0;
    }

    protected override CombinableActionData CombineGuaranteed(CombinableActionData other)
    {
        var behaviourChangeActionData = (BehaviourActionData)other;
        if (Math.Abs(OldValue - behaviourChangeActionData.NewValue) < MathUtils.EPSILON)
            return new BehaviourActionData(behaviourChangeActionData.OldValue, NewValue, Type);

        return new BehaviourActionData(behaviourChangeActionData.NewValue, OldValue, Type);
    }
}