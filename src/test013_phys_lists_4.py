#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gam
from test013_phys_lists_base import create_pl_sim

# create simulation
sim = create_pl_sim()

# keep only ion sources
sim.source_manager.user_info_sources.pop('gamma')

# change physics
p = sim.get_physics_user_info()
p.physics_list_name = 'QGSP_BERT_EMZ'
p.enable_decay = True
mm = gam.g4_units('mm')
cuts = p.production_cuts
cuts.world.gamma = 5 * mm
cuts.world.proton = 1 * mm
cuts.world.electron = -1  # default
cuts.world.positron = 3 * mm
cuts.waterbox.gamma = 2 * mm
cuts.b2.electron = 5 * mm

# initialize
sim.initialize()

print('Phys list cuts:')
print(sim.physics_manager.dump_cuts())

# start simulation
# sim.set_g4_verbose(True)
# sim.apply_g4_command("/tracking/verbose 1")
sim.start()

stats = sim.get_actor('Stats')

# gate_test4_simulation_stats_actor
# Gate mac/main.mac
stats_ref = gam.read_stat_file('./gate/gate_test013_phys_lists/output/stat_4.txt')
is_ok = gam.assert_stats(stats, stats_ref, tolerance=0.1)

gam.test_ok(is_ok)
