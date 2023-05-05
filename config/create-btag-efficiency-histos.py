#!/usr/bin/env python
import itertools
import os
import re
import subprocess
import sys


varial_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Varial')
sys.path.append(varial_dir)
import varial


if len(sys.argv) < 2:
   print "Usage:\n" \
         "    %s <rootfile pattern with full path>\n" \
         "    (all samples are merged)" \
         % os.path.basename(sys.argv[0])
   sys.exit(1)


# definitions
btag_histo_names = (
    "BTagMCEffFlavBPassing",
    "BTagMCEffFlavBTotal",
    "BTagMCEffFlavCPassing",
    "BTagMCEffFlavCTotal",
    "BTagMCEffFlavUDSGPassing",
    "BTagMCEffFlavUDSGTotal",
)


# load and merge histograms
input_pattern = sys.argv[1:]
filter_func = lambda w: w.name in btag_histo_names and not bool(re.search("(TuneCP5|hdamp)(up|down|UP|DOWN)", w.name)) and not bool(re.search("(AToZHToLLTTbar)", w.name))
pipe = varial.gen.open_filter_load(input_pattern, filter_func)
pipe = varial.gen.sort_group_merge(pipe, lambda w: w.name)
pipe = list(pipe)  # pull everything through the generators


# do some checks
if len(pipe) != 6:
    print 'Error. I need histograms for all of these names: ', btag_histo_names
    print 'I only found these: ', list(w.in_file_path for w in pipe)
    exit(-1)

empty_hists = list(w.name for w in pipe if not w.obj.Integral())
if empty_hists:
    print 'Error: integral of histograms "%s" is 0.0! Exit.' % empty_hists
    exit(-1)


# eff histos
pipe = varial.gen.gen_make_eff_graphs(
    pipe,
    postfix_sub='Passing',
    postfix_tot='Total',
    new_postfix='Eff',
    yield_everything=True,
    eff_func=varial.op.div,
)


# write out
outfile = 'BTagMCEfficiencyHists.root'
print 'Saving histograms to file:', outfile
f = varial.ROOT.TFile(outfile, 'RECREATE')
for w in pipe:
    print '  writing histogram:', w.name
    w.obj.SetName(w.name)
    w.obj.Write()
f.Close()
