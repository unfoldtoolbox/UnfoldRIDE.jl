# The Config Struct

When using UnfoldRIDE you will have to supply the function with a `RideConfig` struct. While the Docstring of this is already extensive, we will take the time here to explain every option in detail.

## `sfreq`
The sample frequency of your data. Should be an Integer.

## Component and epoch ranges
`s_range`, `c_range`, and `r_range` indicate the different time-windows around each component, the same way you would choose a time-window to cut your data into epochs.

All component ranges are relative to their respective latency (or estimated latency) in the event structure.

Opposed to that `epoch_range` indicated the entire epoch relative to the S component. This range should encompass all activity from all components.

## `c_estimation_range`
This time window is used for the initial estimation of the C component through peak picking. The range is relative to the S component.

## `tukey_window`
During the cross-correlation a tukey window is applied to each epoch, i.e. every value outside of this range is set to 0. This ensures that a C component can only be found in the range indicated here.

## `formulas`
A Vector containing formulas for the S, R, and C component. These will be used for the Unfold regression model. The formulas will be ignored in the ClassicRIDE algorithm.

## The three heuristics


## Filtering
During the algorithm two filters are applied to the data. By default, both filters are enables (because they are present in the original Matlab implementation), but both can be disabled.

`filtering[1]` enables/disables a 3Hz low-pass filter, which is _only_ applied during the cross-correlation. The filter is applied to deal with the intrinsic oscillatory nature of EEG noise (especially alpha waves).

`filtering[2]` enables/disables a 20Hz low-pass filter which is applied for the (outer) iteration process. This filter is applied to attenuate high frequency noise during latency estimation. The final component estimation uses the original (i.e., unfiltered) data. 