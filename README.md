# Mode-Matching

Aim of this code is to calculate  S matrix of cascaded waveguides with arbitrary width and length for my master thesis.
User only needed to use ```sparam(fmin,fmax,mode_trnc,w,l,s_size)``` function to generate S parameters for a structure.\
Arguments of this function corresponds as folloving:
```
fmin      :minimum frequency
fmax      :maximum frequency
mode_trnc :after which propagating mode to truncate(how many mode to include)
width     :an array that contains width value for each cascaded section
length    :an array that contains length value for each cascaded section
s_size    :number of samples to take in [fmin:fmax] frequency range
```

