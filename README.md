# qaqmc-c
template library implementation for Quasi-Adiabatic Quantum Monte Carlo.

Its implemented for a Nearest Neighbor TFIM based on references:

 1) A. W. Sandvik PHYSICAL REVIEW E 68, 056701 (2003)
 2) C.-W. Liu, A. Polkovnikov, A. W. Sandvik PHYSICAL REVIEW B 87, 174302 (2013)

it implements the 'two-way' evolution.

## Template parameters:
 * int L: linear system size
 * int d: dimension of lattice
 * int M: length of operator string
 * int Nm: Nm <= 2*M+1 number of measurements points during evolution (Nm>0 includes assymetric expectation values).
 * int mstep: number of Monte Carlo sweeps per Monte Carlo step.
 
## Instantiation:
  ```C++
  qaqmc<L,d,M,Nm,mstep> qmc(S_i,S_f,bc,output_file);
  ```

 * `double S_i`: initial 'S' parameter for QA-QMC evolution
 * `double S_f`: final 'S' parameter for QA-QMC evolution
 * `std::string bc`: boundary conditions for operator string, must be one of the following:
   * `"pzpz"` polarized along the positive z-axis on both sides
   * `"pxpx"` polarized along the positive x-axis on both sides
   * `"tr"` tperiodic SSE trace boundary conditions. 
 * `std::string output_file`: string which constains the output file, which the data should be written to. 
 
## public methods of qaqmc:
 * `void write_out(void)`: write the current measurement results to the output file.
 * `void write_out_lock(void)`: writes the current measurement results to output file. the function creates a temporary file which is used to prevent other jobs from writting to thw same file simultaineiously. 
 * `void EQstep(void)`: Perform Monte Carlo step without measurements, typically used when equilibrating simulation.
 * `void MCstep(void)`: Perform Monte Carlo step with measurements, after the completion of this function one should call either `write_out()` or `write_out_lock()`. 
 
