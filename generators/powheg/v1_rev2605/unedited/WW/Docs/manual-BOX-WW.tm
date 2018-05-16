<TeXmacs|1.0.7.3>

<style|article>

<\body>
  <doc-data|<doc-title|The POWHEG-BOX-WW manual>|>

  <section|Introduction>

  The POWHEG-BOX-WW program <cite|noi> can be used to generate the QCD
  production of <math|W<with|math-level|2|<rsup|+>>
  W<rsup|<with|math-level|1|->>> events in hadronic collisions, with the
  <math|W>-bosons decaying into leptons or hadrons, to NLO accuracy in QCD,
  in such a way that matching with a full shower program is possible. In case
  of decays into hadrons, NLO corrections to the decay processes are not
  included. This is unlikely to be necessary: most shower Monte Carlo do
  already a good job in dressing the <math|W> decay with QCD radiation, since
  they have been fit to LEP2 data.

  This generator is based upon the calculation of refs. <cite|Dixon:1998py>,
  <cite|Campbell:1999ah>, <cite|Campbell:2011bn>. The effect of off-shell
  singly resonant graphs is fully included in the calculation. Anomalous
  coupling can also be included.\ 

  \ If the <math|W >-bosons decay into a leptons of the same flavour (e.g.
  <math|e<rsup|+> e<rsup|->>), then the <math|ZZ > production of this signal
  should be considered separately using the POWHEG-BOX-ZZ program.
  Interference between these two process is negligible (see <cite|noi>) and
  is not included. \ This document describes the input parameters that are
  specific to this implementation. The parameters that are common to all
  <with|font-family|tt|POWHEG BOX> implementation are given in the
  <with|font-family|tt|manual-BOX.pdf> document, in the
  <with|font-family|tt|POWHEG-BOX/Docs> directory.

  <section|Generation of events>

  Do<next-line><with|font-family|tt|$ cd POWHEG-BOX/WW><next-line><with|font-family|tt|$
  make pwhg_main><next-line>Then do (for example)<next-line><with|font-family|tt|$
  cd test><next-line><with|font-family|tt|$ ../pwhg_main><next-line>At the
  end of the run, the file <with|font-family|tt|pwgevents.lhe> will contain
  events for <math|W >-pair production in the Les Houches format. In order to
  shower them with <with|font-family|tt|PYTHIA>:<next-line><with|font-family|tt|$
  cd POWHEG-BOX/WW><next-line><with|font-family|tt|$ make
  main-PYTHIA-lhef><next-line><with|font-family|tt|$ cd
  test><next-line><with|font-family|tt|$ ../main-PYTHIA-lhef>

  <section|Input parameters>

  Parameters in <with|font-family|tt|powheg.input> that are specific to
  <math|W W> pair production:<next-line><with|font-family|tt|vdecaymodeWp -11
  \ \ \ ! decay mode of W+ (-1=dbar, -3=sbar, -7=dbar or sbar,<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ !
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -11=e+, -13=mu+,
  -15=tau+.)><next-line><with|font-family|tt|vdecaymodeWm 13 \ \ \ \ ! decay
  mode \ of W- (1=d, 3=s, 7=d or s, <next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ !
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 11=e-,13=mu-,15=tau-.)><with|font-family|tt|<next-line><with|font-family|tt|dronly
  \ \ \ \ \ \ 0 \ \ \ \ \ ! (default 0), if 1 include only double<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! resonant
  contributions><next-line>zerowidth \ \ \ 0 \ \ \ \ \ ! (default 0), if 1
  use on-shell <with|mode|math|W>-bosons only>
  <next-line><with|font-family|tt|runningwidth 0 \ \ \ \ \ ! (default 0), if
  1 use running width> \ \ <next-line><with|font-family|tt|fixedscale \ \ -1
  \ \ \ \ ! (default -1) if \<gtr\>= 0, use fixed scale mu = M_W
  \ \ \ \ \ \ \ \ ><with|font-family|tt|btlscalereal 0 \ \ \ \ \ ! (default
  0) if 1, changes scale for real emission>
  <next-line><with|font-family|tt|btlscalect \ \ 0 \ \ \ \ \ ! (default 0) if
  1, sets scale in counterterm<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ! \ same as Born scale>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ <next-line><with|font-family|tt|withdamp
  \ \ \ \ 1 \ \ \ \ \ ! (default 1) use Born-zero damping factor >

  \;

  The CKM matrix is assumed diagonal. As an
  example<next-line><with|font-family|tt|vdecaymodeWp -7 \ \ \ \ ! decay mode
  of W+ (-1=dbar, -3=sbar, -7=dbar or sbar,<next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ !
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -11=e+, -13=mu+,
  -15=tau+.)><next-line><with|font-family|tt|vdecaymodeWm 13 \ \ \ \ ! decay
  mode \ of W- (1=d, 3=s, 7=d or s, <next-line>
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ !
  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 11=e-,13=mu-,15=tau-.)><next-line>is
  the semileptonic decay <math|\<mu\> <wide|\<nu\>|\<bar\>><rsub|><rsub|\<mu\>>
  + hadrons>.

  \;

  If <with|font-family|tt|zerowidth> is absent or equal to zero, the
  <math|W>-boson's are given finite width. Singly resonant graphs are also
  included by default, unless the <with|font-family|tt|dronly> flag is set to
  1. If zerowidth is set to true, dronly is set to true regardless of what is
  in the powheg.input file. Dynamic widths can be used by setting the
  <with|font-family|tt|runningwidth> flag to 1. As a consequence of the
  unitarity property of the CKM matrix, using PDG values for this matrix
  gives the same result as setting it to a diagonal matrix. For convenience,
  a diagonal CKM matrix is used. Seven anomalous couplings are used:
  <with|font-family|tt|delg1_z>, <with|font-family|tt|de1g1_g>,
  <with|font-family|tt|lambda_z>, <with|font-family|tt|lambda_g>,
  <with|font-family|tt|delk_g>, <with|font-family|tt|delk_z>,
  <with|font-family|tt|tevscale> (see section 3.2 of <cite|noi> and
  references therein for a definition of these). These are set to 0 by
  default, unless a non zero value is given in the
  <strong|><with|font-family|tt|powheg.input> file. Note that in POWHEG,
  these couplings are independent of one another. If there are any relations
  (e.g. due to symmetry considerations) between couplings, these should be
  calculated independently and the input set accordingly.

  If the flag <with|font-family|tt|fixedscale> is set greater than or equal
  to 0, then the factorization and renormalization scales are fixed at
  <with|mode|math|M<rsub|W>> <with|mode|math|>. Otherwise, a dynamic scale of
  the mass of the <math|><with|mode|math|><with|mode|math|W>-boson pair will
  be used. If the flag <with|font-family|tt|btlscalereal> is set to 1, then
  the real contribution will be calculated using the transverse momentum of
  the hardest jet as a scale.

  It is also necessary to employ a Born-zero damping factor, to prevent
  enhancement of the NLO due to small Born terms. For this reason, the flag
  <with|font-family|tt|withdamp> should be set to 1. For further details, see
  appendix B of Ref. <cite|noi>.

  <\bibliography|bib|JHEP|paper.bib>
    <\bib-list|1>
      <bibitem*|1><label|bib-noi>T.<nbsp>Melia, P.<nbsp>Nason,
      R.<nbsp>Rontsch, and G.<nbsp>Zanderighi.

      <bibitem*|2><label|bib-Dixon:1998py>L.<nbsp>J. Dixon, Z.<nbsp>Kunszt,
      and A.<nbsp>Signer, <with|font-shape|italic|Helicity amplitudes for
      O(alpha-s) production of <math|W<rsup|+>W<rsup|->>,
      <math|W<rsup|\<pm\>>Z>, <math|Z*Z>, <math|W<rsup|\<pm\>>\<gamma\>>, or
      <math|Z\<gamma\>> pairs at hadron colliders>,
      <with|font-shape|italic|Nucl.Phys.> <with|font-series|bold|B531> (1998)
      3--23, [<hlink|<with|font-family|tt|hep-ph/9803250>|http://xxx.lanl.gov/abs/hep-ph/9803250>].

      <bibitem*|3><label|bib-Campbell:1999ah>J.<nbsp>M. Campbell and
      R.<nbsp>Ellis, <with|font-shape|italic|An Update on vector boson pair
      production at hadron colliders>, <with|font-shape|italic|Phys.Rev.>
      <with|font-series|bold|D60> (1999) 113006,
      [<hlink|<with|font-family|tt|hep-ph/9905386>|http://xxx.lanl.gov/abs/hep-ph/9905386>].

      <bibitem*|4><label|bib-Campbell:2011bn>J.<nbsp>M. Campbell,
      R.<nbsp>Ellis, and C.<nbsp>Williams, <with|font-shape|italic|Vector
      boson pair production at the LHC>, <hlink|<with|font-family|tt|arXiv:1105.0020>|http://xxx.lanl.gov/abs/arXiv:1105.0020>.
      * Temporary entry *.
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|info-flag|detailed>
    <associate|par-hyphen|normal>
    <associate|preamble|false>
    <associate|sfactor|4>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|7|5>>
    <associate|auto-11|<tuple|6|?>>
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|3|1>>
    <associate|auto-4|<tuple|3|2>>
    <associate|auto-5|<tuple|1|2>>
    <associate|auto-6|<tuple|6|5>>
    <associate|auto-7|<tuple|7|5>>
    <associate|auto-8|<tuple|7|5>>
    <associate|auto-9|<tuple|7|5>>
    <associate|bib-Alioli:2010xd|<tuple|2|?>>
    <associate|bib-Alwall:2006yp|<tuple|1|5>>
    <associate|bib-Boos:2001cv|<tuple|2|6>>
    <associate|bib-Campbell:1999ah|<tuple|3|2>>
    <associate|bib-Campbell:2011bn|<tuple|4|2>>
    <associate|bib-Dixon:1998py|<tuple|2|2>>
    <associate|bib-Frixione:2007nw|<tuple|1|?>>
    <associate|bib-Frixione:2007vw|<tuple|3|?>>
    <associate|bib-Frixione:2007zp|<tuple|3|6>>
    <associate|bib-Mangano:1992jk|<tuple|4|6>>
    <associate|bib-Nason:1988xz|<tuple|5|6>>
    <associate|bib-Nason:1989zy|<tuple|6|6>>
    <associate|bib-Nason:2004rx|<tuple|7|6>>
    <associate|bib-noi|<tuple|1|2>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      noi

      Dixon:1998py

      Campbell:1999ah

      Campbell:2011bn

      noi

      noi

      noi
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Generation
      of events> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Input
      parameters> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>