      t4 = 1/spb(k1,k2)
      t5 = 1/spa(k1,k2)*t4
      t7 = 1/spa(k2,k3)
      t9 = 1/dotp(k1,kq)
      t12 = 1/dotp(k3,kqb)
      t14 = t5*t7*t9*t12
      t15 = spa(k2,q1)
      t16 = spa(k2,q2)
      t17 = t15*t16
      t18 = spa(k2,r1)
      t20 = spb(k1,q1)
      t21 = spb(k1,r2)
      t22 = t20*t21
      t23 = spb(k3,r1)
      t27 = spa(k2,r2)
      t29 = spb(k3,r2)
      t33 = t16**2
      t35 = spb(k1,q2)
      t36 = t35*t21
      t45 = 1/spa(r1,r2)
      t47 = t18**2
      t49 = spa(q1,q2)
      t50 = t20**2
      t51 = t49*t50
      t56 = t9*t12
      t58 = t5*t7*t56*t45
      t59 = t15*t18
      t70 = t16*t18
      t73 = t35*t29
      t78 = 1/spb(k2,k3)
      t80 = m**2
      t83 = 1/(2*t80+2*dotp(kq,kqb))
      t85 = spb(k1,k3)
      t86 = t85**2
      t94 = spb(k3,q1)
      t99 = 1/dotp(k2,k3)
      t100 = t9*t99
      t101 = t5*t100
      t102 = t85*t20
      t104 = t17*t102*t29
      t107 = t33*t85*t73
      t110 = t5*t100*t45
      t120 = t5*t56
      t123 = t4*t83
      t128 = t99*t45
      t132 = t18*t49*t86*t20
      aPMP = t14*t17*t18*t22*t23/4+t14*t17*t27*t22*t29/4+t14*t33*t18*t36
     #*t23/4+t14*t33*t27*t36*t29/4-t14*t45*t15*t47*t51*t23/4-t58*t59*t27
     #*t51*t29/4-t58*t16*t47*t49*t20*t35*t23/4-t58*t70*t27*t49*t20*t73/4
     #+t5*t78*t83*t16*t86*t29-t5*t78*t83*t45*t18*t49*t86*t94+t101*t104/4
     #+t101*t107/4-t110*t59*t49*t102*t94/4-t110*t70*t49*t85*t35*t94/4-t1
     #20*t104/4-t120*t107/4+t123*t99*t16*t86*t21/2-t123*t128*t132/2+t4*t
     #9*t128*t132/4
