fem def para;r;LV_Cubic
fem def coor;r;LV_Cubic
fem def reg;r;LV_Cubic
fem def base;r;LV_Cubic

fem def node;r;/hpc/zwan145/Heart_Failure_Project_STF_Workspace/GeomData/MR_160160/Passive/MR_160160_14
fem def elem;r;/hpc/zwan145/Heart_Failure_Project_STF_Workspace/GeomData/MR_160160/Passive/MR_160160_14

fem def base;r;LV_Cubic_UnitScale
fem update scale_factors unit geometry

fem exp node;Models_UnitSF/DSWall_MR_160160
fem exp elem;Models_UnitSF/DSWall_MR_160160
fem def node;w;Models_UnitSF/DSWall_MR_160160
fem def elem;w;Models_UnitSF/DSWall_MR_160160
fem def elem;w;DSWall_unit

fem quit;
