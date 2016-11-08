fem def para;r;LV_Cubic
fem def coor;r;LV_Cubic
fem def reg;r;LV_Cubic
fem def base;r;LV_Cubic

# Generate surface data points from original model.
fem def node;r;Models_UnitSF/DSWall_STF_01
fem def elem;r;Models_UnitSF/DSWall_STF_01

fem def data;r;Temp
fem def xi;r;Surface_Points_Endo
fem def data;c from_xi
fem def data;w;Models_UnitSF/DSWall_STF_01_SP_Endo
fem exp data;Models_UnitSF/DSWall_STF_01_SP_Endo

fem def data;r;Temp
fem def xi;r;Surface_Points_Epi
fem def data;c from_xi
fem def data;w;Models_UnitSF/DSWall_STF_01_SP_Epi
fem exp data;Models_UnitSF/DSWall_STF_01_SP_Epi

fem reallocate

fem def para;r;LV_Cubic
fem def coor;r;LV_Cubic
fem def reg;r;LV_Cubic
fem def base;r;LV_Cubic

# Project to reconstructed model.
fem def node;r;ReconstructionModels/Reconstructed_STF_01
fem def elem;r;Models_UnitSF/DSWall_STF_01

fem def data;r;Models_UnitSF/DSWall_STF_01_SP_Endo
fem group faces 4,8,12,15,20,24,28,31,36,40,44,47,52,56,60,63 as ENDO
fem def xi;c closest_face ENDO search_start 20
fem list data;ReconstructionModels/STF_01_Endo_Error error
fem export data;ReconstructionModels/STF_01_Endo_Error as EndoError error

fem def data;r;Models_UnitSF/DSWall_STF_01_SP_Epi
fem group faces 5,9,13,16,21,25,29,32,37,41,45,48,53,57,61,64 as EPI
fem def xi;c closest_face EPI search_start 20
fem list data;ReconstructionModels/STF_01_Epi_Error error
fem export data;ReconstructionModels/STF_01_Epi_Error as EpiError error

fem quit
