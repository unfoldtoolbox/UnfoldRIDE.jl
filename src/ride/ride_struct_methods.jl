function response_component_label(cfg::RideConfig)
    @assert length(cfg.component_order) >= 2 "component_order must contain at least anchor and response components"
    return cfg.component_order[2]
end

function component_range(label::Char, cfg::RideConfig)
    return get_component_spec_by_label(label, cfg).range
end

function component_labels(cfg::RideConfig)
    return [comp.label for comp in cfg.components]
end

function get_component_spec_by_label(label::Char, cfg::RideConfig)
    for comp in cfg.components
        if comp.label == label
            return comp
        end
    end
    throw(ArgumentError("No component spec found for label: $label"))
end

function component_events(evts::DataFrame, label::Char)
    return @subset(evts, :event .== label)
end

function anchor_events(evts::DataFrame, cfg::RideConfig)
    return component_events(evts, cfg.anchor_component)
end

function get_component_for_event(raw_event, cfg::RideConfig)
    raw_event_string = string(raw_event)
    for comp in cfg.components
        if any(code -> string(code) == raw_event_string, comp.raw_event_codes)
            return comp
        end
    end
    throw(ArgumentError("No component spec found for raw event code: $raw_event"))
end

function normalize_event_components(evts::DataFrame, cfg::RideConfig)
    @assert :event in names(evts) "Event table must contain an :event column"
    evts_copy = copy(evts)
    if :event_raw ∉ names(evts_copy)
        evts_copy[!, :event_raw] = evts_copy[!, :event]
    end
    evts_copy[!, :event] = [get_component_for_event(ev, cfg).label for ev in evts_copy[!, :event]]
    return evts_copy
end

function adjusted_component_range(comp::ComponentSpec)
    if comp.variable_latency
        return [0, comp.range[2] - comp.range[1]]
    end
    return comp.range
end

function unfold_model_component_entries(cfg::RideConfig; include_variable_latency::Bool = true)
    comps = [comp for comp in cfg.components if include_variable_latency || !comp.variable_latency]
    @assert !isempty(comps) "No component specs available for Unfold model construction"
    return [
        comp.label => (comp.formula, firbasis(adjusted_component_range(comp), cfg.sfreq, ""))
        for comp in comps
    ]
end

function variable_latency_labels(cfg::RideConfig)
    return [comp.label for comp in cfg.components if !isnothing(comp.estimation_range)]
end

function single_variable_latency_label(cfg::RideConfig)
    labels = variable_latency_labels(cfg)
    @assert length(labels) == 1 "Only one variable latency component is supported currently"
    return labels[1]
end