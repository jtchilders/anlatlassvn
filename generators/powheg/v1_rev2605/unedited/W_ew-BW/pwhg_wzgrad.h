        real*8 alpha0,alphas
	common/couplings/alpha0,alphas

	real*8 me,mmu,mtop
	common/masses/me,mmu,mtop

        real*8 mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
	common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs

	real*8 mff,mffs,mii
	common/par_mass/mff,mffs,mii

	real*8 meff
	common/mlight/meff

	real*8 mh,mw,mz,mv,xmh,xmw,xmz
	common/masses_bos/mh,mw,mz,mv,xmh,xmw,xmz

	integer wopt
	real*8 gamv,gamw,gamz
	common/widths/gamv,gamw,gamz,wopt

	real*8 cw,cw2,sw,sw2
	common/mixing/cw,cw2,sw,sw2

	real*8 mu_f,mu_r
	common/scales/mu_f,mu_r

	real*8 deltac,deltas
	common/slicing/deltac,deltas

	integer rep
	real*8 gfermi,alpha
	common/electroweak/gfermi,alpha,rep

	real*8 w2
	common/misc/w2

	integer qcd,mzero,lfc,test(10),qnonr,QED
	common/options/qcd,mzero,lfc,test,qnonr,QED

	integer collcut
	real*8 lambda
	common/cutoffs/lambda,collcut

	complex*16 ieps
	common/cc/ieps

	real*8 flg_wgrad2
	common/flags/flg_wgrad2

	real*8 br_sinv(4,4),rl_sinv(5,5)
	common/invariants/br_sinv,rl_sinv



