# Notebook setup
import os
import sys
import pyrosetta

from pyrosetta import *
from pyrosetta.toolbox import cleaning
from pyrosetta.toolbox.rcsb import pose_from_rcsb
import pyrosetta.rosetta.protocols.surface_docking as sd
from pyrosetta.rosetta.protocols.docking import setup_foldtree
####################
'''
主要目的：
存在一个金属界面，在空间中插入一个蛋白质，不断旋转产生构象，垂直落到距离表面一定距离，再进行微调，
使更多的氨基酸作用到表面上，最后对结构进行一个很小的优化，实现docking

存在的问题：
蛋白与金属的作用没有很好的描述（不知如何设置参数，打分），在最后那一步构象优化始终是硬接触（我更希望是柔性对接），
与MD得到的构象不同。相反，这个方法之前的文章，在磷灰石上的结果很好。

一个比较简单的思路是：
在此基础上，修正最后硬接触部分，使其可以柔性的与表面作用，输出当前结构、n种结构微调的可能性以及score，就很好了

最好的是：提供两种方法，一种快，同上；另一种可以柔性对接，模拟真实过程，但耗时   ————————  与前人文章的意思相同

目标：
得到m种可能的构象接触方式，以及n种可能构象的残基微调（fluctuations），即m*n个可能构象
'''
####################
pyrosetta.init()     #初始化pyrosetta
def num(pose):       #寻找距离界面小于8埃的氨基酸
    total = pose.total_residue()
    attached_AA = 0
    for i in range(2501, total):
        if pose.residue(i).xyz("CA")[2] < 10:
            attached_AA = attached_AA + 1
    return attached_AA

def num1(pose):       #寻找距离界面小于0埃的氨基酸
    total = pose.total_residue()
    attached_AA = 0
    for i in range(2501, total):
        if pose.residue(i).xyz("CA")[2] < 3:
            attached_AA = attached_AA + 1
    return attached_AA
####################

# pose = pose_from_rcsb("1V74")
pose = pose_from_file("gold.pdb")
starting_pose = pose.clone()
cen_pose = pose.clone()
# cen_switch = SwitchResidueTypeSetMover("centroid")    #可将氨基酸结构简化至著有主链
# cen_switch.apply(cen_pose)
starting_cen_pose = cen_pose.clone()

from pyrosetta.rosetta.protocols.docking import setup_foldtree

#添加一个跳跃点？B动A不动，看foldtree
print(pose.fold_tree())
print("before\n")
setup_foldtree(pose,"U_A", Vector1([1]))
print(pose.fold_tree())
print("after\n")
# setup_foldtree(starting_pose,"U_B", Vector1([1]))
# print(starting_pose.fold_tree())
#旋转矩阵，位移向量、中心在哪
jump_num = 1    #很关键
# print(pose.jump(jump_num).get_rotation())
# print('\n')
# print(pose.jump(jump_num).get_translation())

import pyrosetta.rosetta.protocols.rigid as rigid_moves
pert_mover = rigid_moves.RigidBodyPerturbMover(jump_num, 30, 0)

from pyrosetta import PyMOLMover
pymol = PyMOLMover()
pymol.apply(pose)
'''产生足够差异的初始构象'''
pert_mover.apply(pose)
pymol.apply(pose)
pert_mover.apply(pose)
pymol.apply(pose)
pert_mover.apply(pose)
pymol.apply(pose)
pert_mover.apply(pose)
pymol.apply(pose)

'''随机的刚性移动'''
# randomize1 = rigid_moves.RigidBodyRandomizeMover(pose, jump_num, rigid_moves.partner_upstream)
# randomize2 = rigid_moves.RigidBodyRandomizeMover(pose, jump_num, rigid_moves.partner_downstream)
# randomize1.apply(pose)
# pymol.apply(pose)
###########################################
'''针对不同模式的刚性对接'''
# slide = pyrosetta.rosetta.protocols.docking.DockingSlideIntoContact(jump_num)  # for centroid mode
# slide.apply(pose)
slide = pyrosetta.rosetta.protocols.docking.FaDockingSlideIntoContact(jump_num)  # for full-atom mode
slide.apply(pose)
pymol.apply(pose)
# ah = sd.FullatomRelaxMover()
# ah.apply(pose)
movemap = MoveMap()
movemap.set_jump(jump_num, True)
min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
min_mover.movemap(movemap)
scorefxn = get_fa_scorefxn()
min_mover.score_function(scorefxn)
scorefxn(pose)
# print(scorefxn)
min_mover.apply(pose)
# print(pose.jump(jump_num).get_rotation())
# print(pose.jump(jump_num).get_translation())
pymol.apply(pose)
#
'''对蛋白质结构打分，但这个打分不对'''
Prot_trans = rigid_moves.RigidBodyTransMover(pose,jump_num)
sfxn = pyrosetta.rosetta.protocols.loops.get_fa_scorefxn()

pymol.apply(pose)
print("This is the score:", sfxn(pose))
# pymol.send_energy(pose)
###############################################################
while sfxn(pose) > 3560:
    rb_trans = pyrosetta.rosetta.protocols.rigid.RigidBodyTransMover(pose, jump_num, False)
    rb_trans.apply(pose)
    pert_mover.apply(pose)
    pymol.apply(pose)
    slide = pyrosetta.rosetta.protocols.docking.FaDockingSlideIntoContact(jump_num)  # for full-atom mode
    slide.apply(pose)
    sfxn(pose)
    print("This is the score:", sfxn(pose))


# select_heavy_chain = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('U')
# select_light_chain = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('A')
# light_chain = pose.chain('A')
# print(light_chain)


'''对构象进行小的调整，使更多的氨基酸接触蛋白质表面'''
small_mover = rigid_moves.RigidBodyPerturbMover(jump_num, 8, 0)
for i in range(500):
    before_pose = Pose()
    before_pose.assign(pose)
    small_mover.apply(pose)
    slide.apply(pose)
    pymol.apply(pose)
    after_pose = Pose()
    after_pose.assign(pose)
    if num(before_pose) < num(after_pose) and num1(after_pose) < 1:
        before_pose.assign(after_pose)
        pose.assign(after_pose)
    else:
        pose.assign(before_pose)

'''在得到表面最大接触氨基酸后，对接近界面的氨基酸进行优化(相当于柔性对接)'''





'''对非表面的氨基酸能量进行优化'''
pack_score = create_score_function('ref2015')
dock_score = create_score_function('ref2015','docking')
dock_prepack = pyrosetta.rosetta.protocols.docking.DockingPrepackProtocol()
moveable_jump = pyrosetta.rosetta.utility.vector1_int()
moveable_jump.append(1)
dock_prepack.set_partners('U_A')
dock_prepack.set_movable_jumps(moveable_jump)
dock_prepack.set_dock_ppk(True)
dock_prepack.score_and_output("gold.pdb", pose)
dock_prepack.set_sc_min(True)
dock_prepack.set_scorefxn(dock_score)
dock_prepack.set_scorefxn_pack(pack_score)
dock_prepack.apply(pose)
pymol.apply(pose)

##################################################################################################
'''蒙特卡洛对接方式，仍是刚性对接'''
# dock_mcm = pyrosetta.rosetta.protocols.docking.DockMCMProtocol(jump_num, dock_score, pack_score)
# print("here it is")
# dock_mcm.set_first_cycle(50)
# dock_mcm.set_second_cycle(10)
# dock_mcm.set_move_map(movemap)
# dock_mcm.set_task_factory(pyrosetta.rosetta.core.pack.task.TaskFactory())
# dock_mcm.apply(pose)
# pymol.apply(pose)

# scorefxn_high = create_score_function("ref2015.wts", "docking")
# dock_hires = pyrosetta.rosetta.protocols.docking.DockMCMProtocol()
# dock_hires.set_scorefxn(scorefxn_high)
# dock_hires.set_partners("U_A")

# mcm_cycle = pyrosetta.rosetta.protocols.docking.DockMCMCycle()
# mcm_cycle.set_min_type('U_A')
# mcm_cycle.set_move_map(movemap)
# mcm_cycle.set_rot_magnitude(5.0)
# mcm_cycle.set_rtmin(True) # 十分耗时。
# mcm_cycle.set_scmin(True) #
# mcm_cycle.set_task_factory(pyrosetta.rosetta.core.pack.task.TaskFactory())
# pymol.apply(pose)
#################################################################################################
#################################################################################################
'''protein-protein docking'''
# #Low-resolution docking via rosettadock
# scorefxn_low = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("interchain_cen")
# dock_lowres = pyrosetta.rosetta.protocols.docking.DockingLowRes(scorefxn_low, jump_num)
# print(cen_pose.fold_tree())
# setup_foldtree(cen_pose, "A_B", Vector1([1]))
# print(cen_pose.fold_tree())
# dock_lowres.apply(cen_pose)
# print(pyrosetta.rosetta.core.scoring.CA_rmsd(cen_pose, starting_cen_pose))
# print(pyrosetta.rosetta.protocols.docking.calc_Lrmsd(cen_pose, starting_cen_pose, Vector1([1])))
# pymol.keep_history(True)
# pymol.apply(cen_pose)
# pymol.apply(pose)
# 输出pdb
# cen_pose.dump_pdb("cen_pose.pdb")
# pose.dump_pdb("pose.pdb")


'''批处理，输出10个构象'''
# jd = pyrosetta.toolbox.py_jobdistributor.PyJobDistributor("output", 10, scorefxn_low)
# jd.native_pose = starting_cen_pose
#
# while not jd.job_complete:
#     cen_pose.assign(starting_cen_pose)
#     dock_lowres.apply(cen_pose)
#     jd.output_decoy(cen_pose)
#
# pose1 = pose_from_file("output_1.pdb")
# cen_switch = SwitchResidueTypeSetMover("centroid")
# cen_switch.apply(pose1)
# starting_pose1 = pose1.clone()
# setup_foldtree(pose1, "A_B", Vector1([1]))
# print(pose1.fold_tree())
# scorefxn_high = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015.wts", "docking")
# dock_hires = pyrosetta.rosetta.protocols.docking.DockMCMProtocol()
# dock_hires.set_scorefxn(scorefxn_high)
# dock_hires.set_partners("A_B")  # make sure the FoldTree is set up properly
# recover_sidechains = pyrosetta.rosetta.protocols.simple_moves.ReturnSidechainMover(starting_pose1)
# recover_sidechains.apply(pose1)
# print(pyrosetta.rosetta.protocols.surface_docking.SurfaceDockingProtocol().apply(pose))
#################################################################################################################



























































































