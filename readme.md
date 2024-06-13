# Add Permissions for Wham Subprocess

#### MacOS
This command adds executable permissions to the owner, group, and complementary:
~~~
cd path_to_scoresolver
cd wham
cd wham
chmod +x wham
~~~

To verify the command worked, run:
~~~
ls -l wham
~~~

If it says 
~~~
chmod: wham: No such file or directory
~~~

Execute the command to compile wham.c
~~~
make
~~~

Example Output:

~~~
-rwxr-xr-x  1 user  group  12345678 date wham
~~~


# Program
~~~mermaid
---
title: Bind Coefficient Calculation
---
graph LR;
    pdb[("test.pdb")]
    win["windows"]
    mt["metafile.txt & cv data"]
    pmf["pmf.txt"]
    traj["smd_traj.dcd"]
    wl["wham_logs.txt"]
    pdb -- "SMD Pulling Loop" --> traj & win
    win -- "Umbrella Sampling" --> mt
    mt -- "WHAM ANALYSIS" --> pmf & wl
~~~

## Troubleshooting

If the histograms have poor overlap you will need to use more windows and/or reduce the spring constant.

The windows do not need to be linearly spaced and they can have different spring constants.

## Theory

### Definitions

- \( F \): Free energy
- \( k_B \): Boltzmann constant
- \( T \): Temperature
- \( Z \): Canonical partition function
- \( \xi \): Collective Variable (CV)
- \( \mathbf{r} \): Atomic coordinates
- \( P(\xi) \): Probability distribution of the system as a function of CV
- \( \beta \): Inverse temperature (\(\beta = 1 / k_B T\))
- \( U(\mathbf{r}) \): Potential energy of the system
- \( U_i^{\text{bias}} \): Biasing potential for window \( i \)
- \( k \): Force constant for the biasing potential
- \( \xi_i \): Target value of the CV in window \( i \)
- \( P_i^{\text{bias}}(\xi) \): Biased probability distribution for window \( i \)
- \( n_i(\xi) \): Number of samples in bin \( \xi \) for window \( i \)
- \( F_i^{\text{bias}}(\xi) \): Biased free energy in window \( i \)
- \( F(\xi) \): Unbiased free energy (Potential of Mean Force, PMF)

### Equations

1. Free energy in the canonical ensemble:
   \[
   F = -k_B T \ln Z
   \]

2. Probability distribution of the system:
   \[
   P(\xi) = \frac{1}{Z} \int \mathrm{d}\mathbf{r} \, \delta(\xi - \xi(\mathbf{r})) e^{-\beta U(\mathbf{r})}
   \]

3. Free energy as a function of probability:
   \[
   F(\xi) = -k_B T \ln P(\xi) + \text{constant}
   \]

4. Biasing potential for umbrella sampling:
   \[
   U_i^{\text{bias}}(\xi) = \frac{1}{2} k (\xi - \xi_i)^2
   \]

5. Potential energy with biasing potential:
   \[
   U^{\text{total}} = U(\mathbf{r}) + U_i^{\text{bias}}(\xi)
   \]

6. Biased probability distribution:
   \[
   P_i^{\text{bias}}(\xi) = \frac{1}{Z_i} \int \mathrm{d}\mathbf{r} \, \delta(\xi - \xi(\mathbf{r})) e^{-\beta (U(\mathbf{r}) + U_i^{\text{bias}}(\xi))}
   \]

7. Biased free energy:
   \[
   F_i^{\text{bias}}(\xi) = F(\xi) + U_i^{\text{bias}}(\xi)
   \]

8. Combining windows using WHAM:
   \[
   P_i^{\text{bias}}(\xi) = \frac{n_i(\xi)}{\sum_j n_j(\xi) \exp\left(\frac{F_j^{\text{bias}}(\xi) - U_j^{\text{bias}}(\xi)}{k_B T}\right)}
   \]

9. Unbiased free energy profile (PMF):
   \[
   F(\xi) = -k_B T \ln P(\xi)
   \]
