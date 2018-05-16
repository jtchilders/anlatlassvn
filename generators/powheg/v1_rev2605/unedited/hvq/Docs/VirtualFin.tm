<TeXmacs|1.0.7.3>

<style|generic>

<\body>
  <doc-data|<doc-title|Form of the virtual term>>

  The generalization to massive partons is the following. We define
  <math|<with|math-font|cal|I<rsub|l>>> and
  <math|<with|math-font|cal|I><rsub|m>> to be the massless and massive
  sunbsets of <math|<with|math-font|cal|I>>. We also define

  <\equation*>
    \<beta\><rsub|i j>=<sqrt|1-<frac|k<rsub|i><rsup|2>k<rsub|j><rsup|2>|(k<rsub|i>\<cdot\>k<rsub|j>)<rsup|2>>>
    .
  </equation*>

  We have:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<with|math-font|cal|V><rsub|b>>|<cell|=>|<cell|<with|math-font|cal|N><frac|\<alpha\><rsub|s>|2\<pi\>><left|[>-<big|sum><rsub|i\<in\><with|math-font|cal|I><rsub|l>><left|(><frac|1|\<epsilon\><rsup|2>>C<rsub|f<rsub|i>>+<frac|1|\<epsilon\>>\<gamma\><rsub|f<rsub|i>><right|)>
    <with|math-font|cal|B> >>|<row|<cell|>|<cell|>|<cell|-<frac|1|\<epsilon\>><big|sum><rsub|j\<in\><with|math-font|cal|I><rsub|m>>C<rsub|f<rsub|j>>
    <with|math-font|cal|B>>>|<row|<cell|>|<cell|>|<cell|+<frac|2|\<epsilon\>><big|sum><rsub|i,j
    \<in\><with|math-font|cal|I><rsub|l>.
    i\<gtr\>j>log<frac|2k<rsub|i>\<cdot\>k<rsub|j>|Q<rsup|2>>
    \ \ <with|math-font|cal|B><rsub|i j>>>|<row|<cell|>|<cell|>|<cell|+<frac|2|\<epsilon\>><big|sum><rsub|i\<in\>I<rsub|l>,j\<in\><with|math-font|cal|I<rsub|m>>><left|(>log<frac|2k<rsub|i>\<cdot\>k<rsub|j>|Q<rsup|2>>-<frac|1|2>
    log<frac|m<rsub|j><rsup|2>|Q<rsup|2>><right|)>
    <with|math-font|cal|B><rsub|i j>>>|<row|<cell|>|<cell|>|<cell|+<frac|1|\<epsilon\>><big|sum><rsub|i,j
    \<in\><with|math-font|cal|I><rsub|m>, i\<gtr\>j><frac|1|\<beta\><rsub|i
    j>>log<frac|1+\<beta\><rsub|i j>|1-\<beta\><rsub|i j>>
    \ \ <with|math-font|cal|B><rsub|i j> \ \ +
    \ <with|math-font|cal|V><rsub|fin><right|]>,<eq-number><label|eq:Vfin>>>>>
  </eqnarray*>

  where the first term is as in the massless case; the second term arises
  from massive lines emitting and reabsorbing a gluon; the third line is as
  in the massless case. The fourth line is from a light-parton heavy-parton
  interference. The fift line is from heavy-heavy interference.

  These formulae can be inferred from appendix A of the BOX paper. We begin
  by considering the massless case. The divergent part of the soft
  contribution is given in this case by

  <\equation>
    R<rsub|s>=<with|math-font|cal|N><frac|\<alpha\><rsub|s>|2\<pi\>><big|sum><rsub|i<neg|=>j><left|[><frac|1|\<epsilon\><rsup|2>>+<frac|1|\<epsilon\>><left|(>log<frac|Q<rsup|2>|s\<xi\><rsub|c><rsup|2>>-log2
    k<rsub|i>\<cdot\>k<rsub|j>+log 2k<rsub|i><rsup|0>+log
    2k<rsub|j><rsup|0><rsub|><right|)><right|]><with|math-font|cal|B><rsub|i
    j>,
  </equation>

  where we have used eqs. A.31, A.32 and A.33 of the BOX paper (remembering
  also eq. A.12, A.13 and A.14). This yields for the divergent part of the
  soft contribution

  <\equation>
    R<rsub|s>=<with|math-font|cal|N><frac|\<alpha\><rsub|s>|2\<pi\>><left|[><frac|1|\<epsilon\><rsup|2>><big|sum><rsub|i>C<rsub|f<rsub|i>>
    <with|math-font|cal|B>-<frac|1|\<epsilon\>><big|sum><rsub|i<neg|=>j>log<frac|2k<rsub|i>\<cdot\>k<rsub|j>|Q<rsup|2>>
    <with|math-font|cal|B><rsub|i j>+<frac|1|\<epsilon\>><big|sum><rsub|i>C<rsub|f<rsub|i>>
    <left|(>log<frac|Q<rsup|2>|s\<xi\><rsub|c><rsup|2>>+2 log
    <frac|2k<rsub|i><rsup|0>|Q><right|)><with|math-font|cal|B><right|]> .
  </equation>

  These, combined with the soft-collinear divergent remnants, must cancel the
  divergences of the virtual term. Thus, the collinear divergent remnants
  must have the form

  <\equation>
    C<rsub|s>=<with|math-font|cal|N><frac|\<alpha\><rsub|s>|2\<pi\>><left|{><frac|1|\<epsilon\>><big|sum><rsub|i><left|[>-C<rsub|f<rsub|i>>
    <left|(>log<frac|Q<rsup|2>|s\<xi\><rsub|c><rsup|2>>+2 log
    <frac|2k<rsub|i><rsup|0>|E><right|)>+\<gamma\><rsub|f<rsub|i>><right|]>
    <with|math-font|cal|B><right|}>.
  </equation>

  Now we extend this to the massive case. We have

  <\eqnarray*>
    <tformat|<table|<row|<cell|R<rsub|s>>|<cell|=>|<cell|<with|math-font|cal|N><frac|\<alpha\><rsub|s>|2\<pi\>><left|[><big|sum><rsub|i<neg|=>j\<in\><with|math-font|cal|I<rsub|l>>><frac|1|\<epsilon\>><left|(><frac|1|\<epsilon\>>
    -log<frac|2k<rsub|i>\<cdot\>k<rsub|j>|Q<rsup|2>>
    +log<frac|Q<rsup|2>|s\<xi\><rsub|c><rsup|2>>+2 log <frac|2
    k<rsub|i><rsup|0>|Q><right|)><with|math-font|cal|B><rsub|i
    j>>>|<row|<cell|>|<cell|>|<cell|+2<big|sum><rsub|i\<in\><with|math-font|cal|I><rsub|l>,j\<in\><with|math-font|cal|I<rsub|m>>><frac|1|\<epsilon\>><left|(><frac|1|2\<epsilon\>>-
    log<frac|2k<rsub|i>\<cdot\>k<rsub|j>|Q<rsup|2>>+log<frac|Q<rsup|2>|s\<xi\><rsub|c><rsup|2>>+log<frac|2k<rsub|i><rsup|0>|Q>+log<frac|m<rsub|j>|Q><right|)><with|math-font|cal|B><rsub|i
    j>>>|<row|<cell|>|<cell|>|<cell|-<big|sum><rsub|i<neg|=>j\<in\><with|math-font|cal|I<rsub|m>>><frac|1|\<epsilon\>><frac|1|2\<beta\><rsub|i
    j>>log<frac|1+\<beta\><rsub|i j>|1-\<beta\><rsub|i j>>
    \ \ <with|math-font|cal|B><rsub|i j> >>|<row|<cell|>|<cell|>|<cell|+<big|sum><rsub|i\<in\><with|math-font|cal|I<rsub|m>>><frac|1|\<epsilon\>>C<rsub|f<rsub|i>>
    <with|math-font|cal|B><right|]>.<eq-number><label|eq:Rs>>>>>
  </eqnarray*>

  The first line is again from A.32 and A.33. The second line is from A.26,
  A.27. The third line is from A.39 and A.41. The last line is from A.54,
  remembering eq. A.11. The collinear remnants must have the form

  <\equation>
    C<rsub|s>=<with|math-font|cal|N><frac|\<alpha\><rsub|s>|2\<pi\>><left|{><frac|1|\<epsilon\>><big|sum><rsub|i\<in\><with|math-font|cal|I><rsub|l>><left|[>-C<rsub|f<rsub|i>>
    <left|(>log<frac|Q<rsup|2>|s\<xi\><rsub|c><rsup|2>>+2 log
    <frac|2k<rsub|i><rsup|0>|E><right|)>+\<gamma\><rsub|f<rsub|i>><right|]>
    <with|math-font|cal|B><right|}>.
  </equation>

  Summing up (<reference|eq:Rs>) and (<reference|eq:Cs>) we get exactly the
  opposite of the divergent part of eq. (<reference|eq:Vfin>).

  <\equation*>
    \;
  </equation*>

  \;

  <\equation*>
    \;
  </equation*>

  \;

  \;
</body>

<\initial>
  <\collection>
    <associate|sfactor|4>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|eq:Cs|<tuple|6|2>>
    <associate|eq:Rs|<tuple|5|1>>
    <associate|eq:Vfin|<tuple|1|1>>
  </collection>
</references>