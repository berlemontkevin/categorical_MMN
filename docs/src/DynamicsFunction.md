# DynamicsFunction Documentation

!!! info
    Here are all the functions used to perform the numerical integrations of the networks


## Dynamical functions affecting the dendrites

```@docs
dendrite_input_output!
update_dend!
```


## Updating the currents

```@docs
sum_input!
```

## Updating the firing rates

```@docs
current_to_frequency!
update_firing!
update_adaptation!
```

## Updating the synapses

!!! info
    Inhibition onto the dendrites is slower than inhibition elsewhere (Ali and Thomson 2008)

```@docs
synapse_derivative!
current_synapses!
update_syn!
```

## Temporal integration

```@docs
time_step!
full_time_dynamics!
```
