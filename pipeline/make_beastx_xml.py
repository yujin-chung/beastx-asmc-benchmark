#!/usr/bin/env python3
"""
make_beastx_xml.py — emit a BEAST X (BEAST 1.x-format) XML for one alignment.

Model (matched to simulate.py and to the annealedSMC reference):
  - HKY substitution model, kappa FIXED (default 2.0), base freqs FIXED equal.
  - No gamma (optional --gamma-shape for sensitivity).
  - Strict clock, clock.rate FIXED = 1.0  (branch lengths in subs/site;
    the tree is ultrametric in substitution units).
  - Tree prior (--tree-prior):
      yule      Yule (pure-birth) prior, birth rate estimated  [default]
      expbrlen  prior-sensitivity arm: the annealed SMC program's branch-length
                prior, independent Exp(rate) on every branch (default rate 10,
                i.e. "unconstrained:exp(10)"), carried over to the clock tree.
                The product of the 2n-2 exponential densities is
                rate^(2n-2) * exp(-rate * treeLength), so it is implemented as
                an exponential prior with mean 1/rate on the tree length; no
                speciation prior, no birth rate.
Estimated parameters: node heights (+ yule.birthRate under the Yule prior).  This keeps the
substitution model identical to the truth and to aSMC, so topology differences
reflect the exploration strategy, not the model.

Templated on the structure of BEAST X's own examples/Benchmarks/benchmark2.xml
(verified to run on BEAST X v10.5.0 in this project).  fileLog and treeFileLog
are ACTIVE (the smoke test showed the bundled example had them commented out).

Usage:
  python3 make_beastx_xml.py --fasta aln.fasta --out run.xml \
      --chain-length 2000000 --log-every 1000 --file-stem run \
      --kappa 2.0

Then on the server:
  beast -beagle_SSE -threads N -overwrite -seed S run.xml
"""
import argparse, os
from Bio import SeqIO


def load_fasta(path):
    recs = list(SeqIO.parse(path, "fasta"))
    return [(r.id, str(r.seq)) for r in recs]


def xml_escape(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def build_xml(seqs, chain_length, log_every, file_stem, kappa,
              gamma_shape=None, freqs=(0.25, 0.25, 0.25, 0.25),
              tree_prior="yule", brlen_rate=10.0):
    if tree_prior not in ("yule", "expbrlen"):
        raise ValueError("unknown tree prior %r" % tree_prior)
    yule = tree_prior == "yule"
    taxa = [name for name, _ in seqs]
    out = []
    w = out.append
    w('<?xml version="1.0" standalone="yes"?>')
    w('<beast version="10.5">')

    # --- taxa ---
    w('\t<taxa id="taxa">')
    for name in taxa:
        w('\t\t<taxon id="%s"/>' % xml_escape(name))
    w('\t</taxa>')

    # --- alignment ---
    w('\t<alignment id="alignment" dataType="nucleotide">')
    for name, seq in seqs:
        w('\t\t<sequence>')
        w('\t\t\t<taxon idref="%s"/>' % xml_escape(name))
        w('\t\t\t%s' % seq)
        w('\t\t</sequence>')
    w('\t</alignment>')

    # --- patterns ---
    w('\t<patterns id="patterns" from="1" strip="false">')
    w('\t\t<alignment idref="alignment"/>')
    w('\t</patterns>')

    # --- starting tree via constant-size coalescent simulator ---
    w('\t<constantSize id="constant" units="substitutions">')
    w('\t\t<populationSize><parameter id="constant.popSize" value="0.1" lower="0.0"/></populationSize>')
    w('\t</constantSize>')
    w('\t<coalescentSimulator id="startingTree">')
    w('\t\t<taxa idref="taxa"/>')
    w('\t\t<constantSize idref="constant"/>')
    w('\t</coalescentSimulator>')

    # --- tree model ---
    w('\t<treeModel id="treeModel">')
    w('\t\t<coalescentTree idref="startingTree"/>')
    w('\t\t<rootHeight><parameter id="treeModel.rootHeight"/></rootHeight>')
    w('\t\t<nodeHeights internalNodes="true"><parameter id="treeModel.internalNodeHeights"/></nodeHeights>')
    w('\t\t<nodeHeights internalNodes="true" rootNode="true"><parameter id="treeModel.allInternalNodeHeights"/></nodeHeights>')
    w('\t</treeModel>')

    if yule:
        # --- Yule speciation prior ---
        w('\t<yuleModel id="yule" units="substitutions">')
        w('\t\t<birthRate><parameter id="yule.birthRate" value="1.0" lower="0.0"/></birthRate>')
        w('\t</yuleModel>')
        w('\t<speciationLikelihood id="speciation">')
        w('\t\t<model><yuleModel idref="yule"/></model>')
        w('\t\t<speciesTree><treeModel idref="treeModel"/></speciesTree>')
        w('\t</speciationLikelihood>')
    else:
        # --- iid Exp(rate) on every branch == exponential(mean 1/rate) on tree length ---
        w('\t<treeLengthStatistic id="treeLength">')
        w('\t\t<treeModel idref="treeModel"/>')
        w('\t</treeLengthStatistic>')

    # --- strict clock, rate fixed = 1.0 ---
    w('\t<strictClockBranchRates id="branchRates">')
    w('\t\t<rate><parameter id="clock.rate" value="1.0"/></rate>')
    w('\t</strictClockBranchRates>')

    # --- HKY substitution model, kappa fixed, freqs fixed ---
    w('\t<frequencyModel id="freqModel" dataType="nucleotide">')
    w('\t\t<frequencies><parameter id="frequencies" value="%g %g %g %g"/></frequencies>'
      % tuple(freqs))
    w('\t</frequencyModel>')
    w('\t<HKYModel id="hky">')
    w('\t\t<kappa><parameter id="kappa" value="%g"/></kappa>' % kappa)
    w('\t\t<frequencies><frequencyModel idref="freqModel"/></frequencies>')
    w('\t</HKYModel>')

    # --- site model ---
    w('\t<siteModel id="siteModel">')
    w('\t\t<substitutionModel><HKYModel idref="hky"/></substitutionModel>')
    if gamma_shape is not None:
        w('\t\t<gammaShape gammaCategories="4"><parameter id="alpha" value="%g" lower="0.0"/></gammaShape>'
          % gamma_shape)
    w('\t</siteModel>')

    # --- tree likelihood (BEAGLE) ---
    w('\t<treeLikelihood id="treeLikelihood" useAmbiguities="false">')
    w('\t\t<patterns idref="patterns"/>')
    w('\t\t<treeModel idref="treeModel"/>')
    w('\t\t<siteModel idref="siteModel"/>')
    w('\t\t<branchRateModel><strictClockBranchRates idref="branchRates"/></branchRateModel>')
    w('\t</treeLikelihood>')

    # --- operators ---
    w('\t<operators id="operators">')
    # topology moves
    w('\t\t<subtreeSlide size="0.05" gaussian="true" weight="15"><treeModel idref="treeModel"/></subtreeSlide>')
    w('\t\t<narrowExchange weight="15"><treeModel idref="treeModel"/></narrowExchange>')
    w('\t\t<wideExchange weight="3"><treeModel idref="treeModel"/></wideExchange>')
    w('\t\t<wilsonBalding weight="3"><treeModel idref="treeModel"/></wilsonBalding>')
    # node heights
    w('\t\t<uniformOperator weight="30"><parameter idref="treeModel.internalNodeHeights"/></uniformOperator>')
    w('\t\t<scaleOperator scaleFactor="0.75" weight="3"><parameter idref="treeModel.rootHeight"/></scaleOperator>')
    # birth rate
    if yule:
        w('\t\t<scaleOperator scaleFactor="0.75" weight="3"><parameter idref="yule.birthRate"/></scaleOperator>')
    if gamma_shape is not None:
        w('\t\t<scaleOperator scaleFactor="0.75" weight="1"><parameter idref="alpha"/></scaleOperator>')
    w('\t</operators>')

    # --- MCMC ---
    w('\t<mcmc id="mcmc" chainLength="%d" autoOptimize="true" operatorAnalysis="true">'
      % chain_length)
    w('\t\t<joint id="joint">')
    w('\t\t\t<prior id="prior">')
    if yule:
        w('\t\t\t\t<speciationLikelihood idref="speciation"/>')
        w('\t\t\t\t<oneOnXPrior><parameter idref="yule.birthRate"/></oneOnXPrior>')
    else:
        w('\t\t\t\t<exponentialPrior mean="%.10g" offset="0.0"><treeLengthStatistic idref="treeLength"/></exponentialPrior>'
          % (1.0 / brlen_rate))
    w('\t\t\t</prior>')
    w('\t\t\t<likelihood id="likelihood">')
    w('\t\t\t\t<treeLikelihood idref="treeLikelihood"/>')
    w('\t\t\t</likelihood>')
    w('\t\t</joint>')
    w('\t\t<operators idref="operators"/>')
    # screen log
    w('\t\t<log id="screenLog" logEvery="%d">' % (log_every * 10))
    w('\t\t\t<column label="Joint" dp="4" width="12"><joint idref="joint"/></column>')
    w('\t\t\t<column label="Likelihood" dp="4" width="12"><likelihood idref="likelihood"/></column>')
    w('\t\t\t<column label="rootHeight" dp="4" width="12"><parameter idref="treeModel.rootHeight"/></column>')
    w('\t\t</log>')
    # file log
    w('\t\t<log id="fileLog" logEvery="%d" fileName="%s.log" overwrite="true">'
      % (log_every, file_stem))
    w('\t\t\t<joint idref="joint"/>')
    w('\t\t\t<prior idref="prior"/>')
    w('\t\t\t<likelihood idref="likelihood"/>')
    w('\t\t\t<parameter idref="treeModel.rootHeight"/>')
    if yule:
        w('\t\t\t<parameter idref="yule.birthRate"/>')
    else:
        w('\t\t\t<treeLengthStatistic idref="treeLength"/>')
    w('\t\t\t<treeLikelihood idref="treeLikelihood"/>')
    w('\t\t</log>')
    # tree log
    w('\t\t<logTree id="treeFileLog" logEvery="%d" nexusFormat="true" fileName="%s.trees" sortTranslationTable="true">'
      % (log_every, file_stem))
    w('\t\t\t<treeModel idref="treeModel"/>')
    w('\t\t\t<joint idref="joint"/>')
    w('\t\t</logTree>')
    w('\t</mcmc>')

    w('\t<report>')
    w('\t\t<property name="timer"><mcmc idref="mcmc"/></property>')
    w('\t</report>')
    w('</beast>')
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--file-stem", default=None,
                    help="stem for .log/.trees (default: out basename)")
    ap.add_argument("--chain-length", type=int, default=2000000)
    ap.add_argument("--log-every", type=int, default=1000)
    ap.add_argument("--kappa", type=float, default=2.0)
    ap.add_argument("--gamma-shape", type=float, default=None)
    ap.add_argument("--tree-prior", choices=("yule", "expbrlen"), default="yule")
    ap.add_argument("--brlen-rate", type=float, default=10.0,
                    help="rate of the exponential branch-length prior (expbrlen only)")
    args = ap.parse_args()

    seqs = load_fasta(args.fasta)
    stem = args.file_stem or os.path.splitext(os.path.basename(args.out))[0]
    xml = build_xml(seqs, args.chain_length, args.log_every, stem,
                    args.kappa, gamma_shape=args.gamma_shape,
                    tree_prior=args.tree_prior, brlen_rate=args.brlen_rate)
    open(args.out, "w").write(xml)
    print("wrote %s  (ntax=%d, chainLength=%d, stem=%s)"
          % (args.out, len(seqs), args.chain_length, stem))


if __name__ == "__main__":
    main()
