# Koopmans to Yambo (k2y)

Convert your Koopmans eigenvalues data into yambo-readable ndb.QP databases.
This tool is mainly used to interface `kcw.x` results and `yambo`, e.g. for BSE calculations.
This allows to skip heavy GW calculations.
It can also be used to print the k-points needed by yambo to compute koopmans eigenvalues via interpolation (effectively skipping any additional Koopmans run).

The `examples` folder contains some usage examples for the code:

- `01_generate_ndbQP.py`: generate the ndb.QP starting from ns.db1, a template ndb.QP and koopmans output file;
- `02_generate_ndbQP_aiida.py`: generate the ndb.QP `SinglefileData` starting from a `YamboCalculation` and a `KcwCalculation`;
- `03_produce_kpoints.py`: print the list of kpoints to be interpolated in koopmans;