
       integer maxpa,nwpawc
       parameter (nwpawc=10000000)
       parameter (maxpa=900)
       integer h(nwpawc)

       real E_gn(maxpa),Px_gn(maxpa),Py_gn(maxpa),Pz_gn(maxpa),
     +      Vt_gn(maxpa),Vx_gn(maxpa),Vy_gn(maxpa),Vz_gn(maxpa),
     +	    m_gn(maxpa)

       integer id_gn(maxpa),ngen,istd_gn(maxpa),jda1_gn(maxpa),
     +	       jda2_gn(maxpa),jmo1_gn(maxpa),jmo2_gn(maxpa)
	

       common /pawc/ h
       common /particle/ ngen,E_gn,Px_gn,Py_gn,Pz_gn,Vt_gn,Vx_gn,
     +                   Vy_gn,Vz_gn,id_gn,istd_gn,jda1_gn,
     +                   jda2_gn,jmo1_gn,jmo2_gn,m_gn
     			 

