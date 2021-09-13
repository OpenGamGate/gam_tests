#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gam
import gam_g4
import uproot

# create the simulation
sim = gam.Simulation()

# main options
ui = sim.user_info
ui.g4_verbose = False
ui.visu = False
ui.number_of_threads = 1  # FIXME

# units
m = gam.g4_units('m')
cm = gam.g4_units('cm')
keV = gam.g4_units('keV')
mm = gam.g4_units('mm')
Bq = gam.g4_units('Bq')

# world size
world = sim.world
world.size = [2 * m, 2 * m, 2 * m]

# material
sim.add_material_database('./data/GateMaterials.db')

# fake spect head
waterbox = sim.add_volume('Box', 'SPECThead')
waterbox.size = [55 * cm, 42 * cm, 18 * cm]
waterbox.material = 'G4_AIR'

# crystal
crystal1 = sim.add_volume('Box', 'crystal1')
crystal1.mother = 'SPECThead'
crystal1.size = [0.5 * cm, 0.5 * cm, 2 * cm]
crystal1.translation = None
crystal1.rotation = None
crystal1.material = 'NaITl'
start = [-25 * cm, -20 * cm, 4 * cm]
size = [100, 40, 1]
# size = [100, 80, 1]
tr = [0.5 * cm, 0.5 * cm, 0]
crystal1.repeat = gam.repeat_array('crystal1', start, size, tr)
crystal1.color = [1, 1, 0, 1]

# additional volume
crystal2 = sim.add_volume('Box', 'crystal2')
crystal2.mother = 'SPECThead'
crystal2.size = [0.5 * cm, 0.5 * cm, 2 * cm]
crystal2.translation = None
crystal2.rotation = None
crystal2.material = 'NaITl'
start = [-25 * cm, 0 * cm, 4 * cm]
size = [100, 40, 1]
tr = [0.5 * cm, 0.5 * cm, 0]
crystal2.repeat = gam.repeat_array('crystal2', start, size, tr)
crystal2.color = [0, 1, 0, 1]

# physic list
p = sim.get_physics_user_info()
p.physics_list_name = 'G4EmStandardPhysics_option4'
p.enable_decay = False
cuts = p.production_cuts
cuts.world.gamma = 0.01 * mm
cuts.world.electron = 0.01 * mm
cuts.world.positron = 1 * mm
cuts.world.proton = 1 * mm

# default source for tests
source = sim.add_source('Generic', 'Default')
source.particle = 'gamma'
source.energy.mono = 140.5 * keV
source.position.type = 'sphere'
source.position.radius = 4 * cm
source.position.translation = [0, 0, -15 * cm]
source.direction.type = 'momentum'
source.direction.momentum = [0, 0, 1]
source.activity = 2000 * Bq / ui.number_of_threads

# add stat actor
sim.add_actor('SimulationStatisticsActor', 'Stats')

# hits collection
hc = sim.add_actor('HitsCollectionActor', 'hc')
# hc.mother = [crystal1.name, crystal2.name]  # FIXME
hc.mother = [crystal1.name, crystal2.name]
hc.output = 'hits.root'
hc.branches = ['KineticEnergy', 'PostPosition', 'TotalEnergyDeposit', 'GlobalTime', 'VolumeName']
# hc.branches = ['KineticEnergy', 'PostPosition', 'TotalEnergyDeposit', 'GlobalTime']

# create G4 objects
sim.initialize()

# start simulation
sec = gam.g4_units('second')
sim.run_timing_intervals = [[0, 1 * sec]]
sim.start()

# stat
stats = sim.get_actor('Stats')
print(stats)
stats_ref = gam.read_stat_file('./gate/gate_test025_hits_collection/output/stat.txt')
is_ok = gam.assert_stats(stats, stats_ref, tolerance=0.05)

# root
ref_hits = uproot.open('./gate/gate_test025_hits_collection/output/hits.root')['Hits']
rn = ref_hits.num_entries
ref_hits = ref_hits.arrays(library="numpy")
print(rn, ref_hits.keys())

hits = uproot.open('hits.root')['Hits']
n = hits.num_entries
hits = hits.arrays(library="numpy")
print(n, hits.keys())

diff = gam.rel_diff(float(rn), n)
gam.print_test(diff < 5, f'Nb values: {rn} {n} {diff:.2f}%')

hc = sim.get_actor('hc')
ab = gam_g4.GamBranch.GetAvailableBranches()
for b in ab:
    print(b.fBranchName, b.fBranchType)

tree = hc.GetHits()
print('Tree', tree.fTreeName)
for b in tree.fBranches:
    print('Branch', b.fBranchName, b.size())

# how to compare hits ?
# 1) branch names correspondence, nb value % diff
# 2) if float, compare mean max min std ?
for k in ref_hits:
    is_ok = gam.assert_tree_branch(ref_hits[k], k, hits)
