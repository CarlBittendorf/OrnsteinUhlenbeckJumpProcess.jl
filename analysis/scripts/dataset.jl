function (; var"data#BPD Affect & Self-Esteem Forms")
    @chain var"data#BPD Affect & Self-Esteem Forms" begin
        transform(
            :Participant => ByRow(x -> parse(Int, string(x)[1:6])),
            [:FormTrigger, :FormStart] =>
                ByRow((t, s) -> ismissing(t) ? s : t) => :FormTrigger,
            [:MDMQContentMoment, :MDMQUnwellMoment] =>
                ByRow((c, u) -> (abs(c - 100) + u) / 200) => :Valence,
            :EventValence => ByRow(x -> (x - 50) / 5000);
            renamecols = false
        )

        groupby([:Participant, :Group])
        transform(:FormTrigger => enumerate_days => :Day)

        subset(:Day => ByRow(x -> x <= 4))

        dropmissing(:EventValence)
        transform(
            :EventValence => ByRow(x -> x > 0 ? x : 0.0) => :PositiveEventIntensity,
            :EventValence => ByRow(x -> x < 0 ? abs(x) : 0.0) => :NegativeEventIntensity
        )

        select(:Participant, :Group, :Platform, :FormStart, :Valence,
            :PositiveEventIntensity, :NegativeEventIntensity)
        dropmissing

        groupby([:Participant, :Group])
        subset(
            # use only participants with at least 36 answered queries (75% compliance)...
            :FormStart => (x -> length(x) >= 36),

            # ...variance greater than 0...
            :Valence => (x -> var(x) > 0),

            # ...and with at least two positive and two negative events
            [:PositiveEventIntensity, :NegativeEventIntensity] .=>
                (x -> count(!isequal(0), x) >= 2);
            ungroup = false
        )
        transform(:FormStart => hours_elapsed => :HoursSinceObservation)

        transform(:HoursSinceObservation => ByRow(x -> x / 2) => :HoursSinceEvent)
    end
end