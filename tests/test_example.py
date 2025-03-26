from CoREMOF.calculation import Zeopp


example_cif = "./IRMOF-1.cif"

results_strinfo = Zeopp.FrameworkDim(structure=example_cif,
                                        high_accuracy=True,
                                        prefix="test_strinfo") 
print(results_strinfo["Dimention"])
   