
from k2y.k2y import KcwQpDatabaseGenerator

spin = True  # set to True if the KCW calculation is spin-polarized

aiida_node = 77525

if not aiida_node:
    converter = KcwQpDatabaseGenerator(
        ns_db1="/path/to/ns.db1",
        #template_QP_path="/path/to/template.QP"
        spin=spin
        )
else:
    converter = KcwQpDatabaseGenerator.from_aiida(
        yambo_node_pk=aiida_node,  # Replace with actual pk
        kcw_node_pk=None,    # Replace with actual pk or None
        #template_QP_path="/path/to/template.QP"
        spin=spin
        )
    
converter.set_koopmans_eval(path="kc.kho") # not needed if you are using AiiDA

# we need kpoints from the pw input
converter.set_kpoints_from_pwinput("../1_wannier/pwnscf.in")

"""
if you need to have a smaller subset of bands, to be compatible with Wannier90 eig file:

converter.eigenvalues_KS = converter.eigenvalues_KS[:,:252] # where 252 is num_bands of wann, i.e. of the eig file.
converter.eigenvalues_KI = converter.eigenvalues_KI[:,:252]

"""

converter.generate_mappings()

converter.verify_mappings(k_index=1, top_valence=10) # adjust k_index and top_valence as needed

converter.generate_QP_db("out.QP")

if aiida_node: 
    # if you are using AiiDA, generate the SinglefileData
    # you can also use the converter.generate_QP_db_SinglefileData() method
    new_db = converter.generate_SinglefileData_from_file("out.QP")
    new_db.store()
    print(f"stored {new_db}")


#76693