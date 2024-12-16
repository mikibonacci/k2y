from typing import Union
import netCDF4 as nc
from ase.units import Ha
from ase import io
import xarray
import itertools

#from importlib.path import Path

class KcwQpDatabaseGenerator:
    """
    class to manipulate ns.db1, template.QP and KI eigenvalues to obtain an kcw.QP file.
    we have also a pre-processing method to generate the list of k-points needed in the kcw interpolation,
    reading the ns.db1.
    
    TODO use pydantic model and understand which should be a staticmethod, classmethod...
    """
    def __init__(self,ns_db1:str=None, template_QP_path:str=None):
        """
        ns_db1: path to the ns.db1 netcdf yambo file. This is used to generate the ndb_KI.QP and the K-points list for the interpolation of the
        KI eigenvalues.
        template_QP: path to the template.QP db used to generate the ndb_KI.QP.
        """
        
        if ns_db1:
            print("reading ns.db1")
            self.ns_db1 = xarray.open_dataset(ns_db1,engine='netcdf4')
            
        if template_QP_path:
            print("storing the path of template.QP")
            self.template_QP_path = template_QP_path
    
     
    def produce_kpoints_for_interpolation(self, filename=None):
        """Writes the Kpoints card for kcw.x interpolation. Kpoints are in cartesian coordinates in unit of 2pi/alat.
        TBO (to be optimized).

        Args:
            ns (xarray.core.dataset.Dataset): ns.db1 database of the KS dense mesh
            filename (str, optional): provide the filename path if you want also to write the resulting card, not just print it. Defaults to None.
        """
        
        ns = self.ns_db1
        
        print("K_POINTS tpiba_b")
        nk = len(ns["K-POINTS"].values[0])
        print(f"{nk}")
        for k in range(nk):
            print(f'{ns["K-POINTS"].values[0,k]*2}  {ns["K-POINTS"].values[1,k]*2}  {ns["K-POINTS"].values[2,k]*2} 0')
            
        if filename:
            with open(filename,'w') as file:
                file.write("K_POINTS tpiba_b")
                file.write(f"{nk}\n")
                for k in range(nk):
                    file.write(f'{ns["K-POINTS"].values[0,k]*2}  {ns["K-POINTS"].values[1,k]*2}  {ns["K-POINTS"].values[2,k]*2} 0')
                    
        return 
    
    def set_koopmans_eval(self, path:str=None):
        
        output = io.read(path)
        self.eigenvalues_KI = np.array(output.calc.results["eigenvalues"])
        
        return
    
    def generate_mappings(self,):
        
        ns = self.ns_db1
        self.eigenvalues_KS = ns.variables["EIGENVALUES"].values[0] # not exactly, but we can map.
        
        # use the KI:
        eigenvalues = self.eigenvalues_KI

        print(f"shape eigenvalues: {np.shape(eigenvalues)}")
        
        self.kpoints = ns.variables["K-POINTS"].values[:,:np.shape(eigenvalues)[0]] # exactly as in ndb.QP


        bands = [1,np.shape(eigenvalues)[0]*np.shape(eigenvalues)[1]]
        kpoints_ind = [1,np.shape(self.kpoints)[1]]
        reshaped_eval = eigenvalues.reshape(bands[-1])/Ha
        reshaped_eval_KS = self.eigenvalues_KS[:np.shape(eigenvalues)[0],:np.shape(eigenvalues)[1]].reshape(bands[-1])
        QP_kpts = self.kpoints # [[x],[y],[z]]
        QP_Eo = reshaped_eval_KS # [E_KS]
        QP_Z = [[1,0]]*(bands[1]-bands[0]+1) # [[Z.real,Z.imag]] <--- we don't care about this.
        table_bands = [b_ind for b_ind in range(1,np.shape(eigenvalues)[1]+1)]*np.shape(self.kpoints)[1]
        

        nested_list = [[k]*np.shape(eigenvalues)[1] for k in range(1,1+np.shape(self.kpoints)[1])]
        table_kpoints = list(itertools.chain(*nested_list))

        QP_table = [table_bands,table_bands,table_kpoints]
        
        self.table = QP_table
        
        QP_E = [[E,0] for E in reshaped_eval] # [[E.real,E.imag]]
        PARS = np.ma.array(
            np.array([np.shape(eigenvalues)[1],np.shape(eigenvalues)[0],np.shape(eigenvalues)[0]*np.shape(eigenvalues)[1],26,-1,-1]),
            mask = [False, False, False, False,  True,  True])

        self.mapped_vars = {
        'QP_QP_@_state_1_b_range': [1,np.shape(eigenvalues)[1]],
        'QP_QP_@_state_1_K_range': [1,kpoints_ind[-1]],
        "QP_kpts": QP_kpts,
        "QP_E": QP_E,
        "QP_Eo": QP_Eo,
        "QP_Z": QP_Z,
        "QP_table": QP_table,
        "PARS":PARS,
        }

        # This mapping is done onto template.QP. You should always use it!!!
        self.mapped_dims = {
            'D_0000000064': [f"D_{str(len(QP_Eo)).zfill(10)}",len(QP_Eo)], # 
            'D_0000000008': [f"D_{str(len(QP_kpts[0])).zfill(10)}",len(QP_kpts[0])], # 
            #'D_0000000100': [f"D_{str(mapped_vars['QP_QP_@_state_1_b_range'][1]).zfill(10)}",mapped_vars['QP_QP_@_state_1_b_range'][1]], # 
        }
        
        print("mapping successfully generated.")
        
        return
    
    def generate_QP_db(self,output_filename:str="output.QP"):
        
        print(f"Generating Koopmans quasiparticle database...")
        
        # 1 Open the template file in read-only mode
        template_nc = nc.Dataset(self.template_QP_path,'r')
        new_db = nc.Dataset(output_filename,'w', format='NETCDF4')
        
        # 2.1 Copy dimesions
        for dim_name in template_nc.dimensions:
            if dim_name not in new_db.dimensions:
                if dim_name in self.mapped_dims.keys():
                    new_db.createDimension(self.mapped_dims[dim_name][0], self.mapped_dims[dim_name][1])
                else:
                    new_db.createDimension(dim_name, len(template_nc.dimensions[dim_name]))

        # 2.2 Copy variables
        for var_name in template_nc.variables:
            #print(var_name)
            if var_name not in new_db.variables:
                var = template_nc.variables[var_name]
                new_dims = []
                for v in var.dimensions:
                    if v in self.mapped_dims.keys():
                        new_dims.append(self.mapped_dims[v][0])
                    else:
                        new_dims.append(v)
                new_db.createVariable(var_name, var.dtype, new_dims)
                new_db.variables[var_name][:] = self.mapped_vars[var_name] if var_name in self.mapped_vars.keys() else var[:]

        # 3 Copy global attributes
        for attr_name in template_nc.ncattrs():
            new_db.setncattr(attr_name, template_nc.getncattr(attr_name))
        
        new_db.close()
        template_nc.close()
        
        print("...done.")
        
        return
    
"""Usage example:

converter = KcwQpDatabaseGenerator(
    ns_db1="/path/to/ns.db1",
    template_QP_path="/path/to/template.QP"
    )
    
converter.set_koopmans_eval(path="/path/to/kc.kho")

converter.generate_mappings()

converter.generate_QP_db("out.QP")

#converter.produce_kpoints_for_interpolation() # to produce the k-points card for the interpolation, in the kcw.x run...

"""