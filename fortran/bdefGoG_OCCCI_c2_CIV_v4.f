      program bloomdefGoG_OCCCI_allcase_CIV_suite_c2
CCC---------------------------------------------------------------------
CCC
CC declarations of variables
CC =========================
C
      IMPLICIT NONE
CC
CC    jpi = longitude
CC    jpj = latitude
CC    jpk = depth
CC    jit = number of time step
CC
      INTEGER jpi,jpj,jpk,jit,nvar,nyr
      PARAMETER(jpi=481,jpj=361,jpk=1,jit=828,nvar=4)
      PARAMETER(nyr=17)
      INTEGER ji,jj,jk,it,ii,iit,nv,nn
CC
      REAL wflx(jpi,jpj,jpk,jit),med5(jpi,jpj,jpk,jpk)
      INTEGER dmax(jpi,jpj,jpk,nyr)
CC
      REAL lat(jpj), lon(jpi), depth(jpk), time(jit)
      INTEGER fix
      REAL fix2
C
CC      
      INCLUDE 'netcdf.inc'
c
      INTEGER ncid, iret, varid, ncod, vorid(nvar)
      INTEGER ncid2, varid2, ncid3, varid3
      INTEGER varlat, varlon, vardep, vartim
      INTEGER varlat2, varlon2, vardep2, vartim2
      INTEGER vorlat, vorlon, vordep, vortim
      INTEGER nvdims
      PARAMETER(nvdims = 4)
      INTEGER vardims(nvdims),start(nvdims),count(nvdims),stride(nvdims)
      INTEGER status

CC
      CHARACTER*30 filein(nyr)
cc
      DATA start/1,1,1,1/
      DATA count/jpi,jpj,jpk,jit/
      DATA stride/1,1,1,1/
      DATA fix,fix2/-1000,-1.e+34/
CC
      DATA filein/"bddef98.nc","bddef99.nc","bddef00.nc",
     +     "bddef01.nc","bddef02.nc","bddef03.nc","bddef04.nc",
     +     "bddef05.nc","bddef06.nc","bddef07.nc","bddef08.nc",
     +     "bddef09.nc","bddef10.nc","bddef11.nc","bddef12.nc",
     +     "bddef13.nc","bddef14.nc"/
CC
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  new variables for sorting
c
      INTEGER last, kt, ptr, first, nit, counter, nber
      INTEGER ngd(jpi,jpj), nbd(jpi,jpj)
      INTEGER date_start(jpi,jpj), date_end(jpi,jpj)
      INTEGER datemaxori(jpi,jpj), duration(jpi,jpj) 
      INTEGER date_max(jpi,jpj)
      INTEGER datestartori(jpi,jpj),dateendori(jpi,jpj)
      REAL med0(jpi,jpj,jit), thres(jpi,jpj)
      REAL wflxsort(jpi,jpj,jpk,jit) 
      REAL maxval(jpi,jpj)
      DATA nit/46/
c
c For each year
      INTEGER tstart(nyr),jits(nyr),jan(nyr)
      
      DATA jan/46,92,138,184,230,276,322,368,414,460,506,552,598,644,
     +     690,736,782/
c
      DATA tstart/47,93,139,185,231,277,323,369,415,461,507,553,599,
     +     645,691,737,773/
      DATA jits/92,138,184,230,276,322,368,414,460,506,552,598,644,
     +     690,736,772,828/
      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc loop over years
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       DO nn = 1, nyr

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  open the new netCDF file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
       write(*,*) "opening filein ",filein(nn)
       status = nf_create(filein(nn),nf_clobber,ncod)
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define attributes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_def_dim(ncod,'LONGITUDE',jpi,vardims(1))
       status = nf_def_dim(ncod,'LATITUDE',jpj,vardims(2))
       status = nf_def_dim(ncod,'DEPTH',jpk,vardims(3))
       status = nf_def_dim(ncod,'TIME',1,vardims(4)) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_def_var(ncod,'date_start',nf_int,nvdims,
     +      vardims,vorid(1))
       status = nf_def_var(ncod,'date_max',nf_int,nvdims,
     +      vardims,vorid(2))
       status = nf_def_var(ncod,'date_end',nf_int,nvdims,
     +      vardims,vorid(3))
       status = nf_def_var(ncod,'duration',nf_int,nvdims,
     +      vardims,vorid(4))
C                   
       status=nf_def_var(ncod,'LONGITUDE',nf_float,1,vardims(1),vorlon)
       status=nf_def_var(ncod,'LATITUDE',nf_float,1,vardims(2),vorlat)
       status=nf_def_var(ncod,'DEPTH',nf_float,1,vardims(3),vordep)
       status=nf_def_var(ncod,'TIME',nf_float,1,vardims(4),vortim)
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  create attributes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C  for output
C
       status = nf_put_att_text(ncod, vorlon, 'units', 6, 'degrees')
       status = nf_put_att_text(ncod, vorlat, 'units', 6, 'degrees')
       status = nf_put_att_text(ncod, vordep, 'units', 6, 'meters')
C
C for ocean data:
C
       status = nf_put_att_text(ncod, vordep, 'positive', 4, 'down')
       status = nf_put_att_text(ncod, vortim, 'units', 4, 'year')
C      
C
       DO nv = 1, nvar
       status = nf_put_att_text(ncod,vorid(nv),'units',11,'weeks')
       status = nf_put_att_int(ncod, vorid(nv),'_FillValue', nf_int,
     +      1,fix)
       ENDDO
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  end of define mode
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_enddef(ncod)
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  open the netCDF input file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C NOTE: the variable status should always return a 0 value. When this
C is not 0, then the operation failed (check the name of the file and
C variables, watch for capital letters)
C
      status = nf_open(
     +      "OCCCI_v3-8DAY-97_14_9km_pheno_fav8_sbx3weeks.nc",
     +      nf_nowrite,ncid)
      status = nf_open("mergedbdef_GoG_OCCCI_case_CIV.nc",
     +      nf_nowrite,ncid2)
      status = nf_open("bmed5_GoG_OCCCI_8d.nc",
     +      nf_nowrite,ncid3)

      write(*,*) "File open", status
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  inquire attributes for the grid and variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_inq_varid(ncid,'LONGITUDE',varlon)
       if (status.ne.0) write(*,*) "Problem in varlon",status
C
       status = nf_inq_varid(ncid,'LATITUDE',varlat)
       if (status.ne.0) write(*,*) "Problem in varlat",status
C
       status = nf_inq_varid(ncid,'DEPTH',vardep)
       if (status.ne.0) write(*,*) "Problem in vardep",status
C
       status = nf_inq_varid(ncid,'TIME',vartim)
       if (status.ne.0) write(*,*) "Problem in vartim",status
C
       status = nf_inq_varid(ncid,'CHL',varid)
       if (status.ne.0) write(*,*) "Problem in varid",status
       status = nf_inq_varid(ncid2,'date_max',varid2)
       if (status.ne.0) write(*,*) "Problem in varid2",status
       status = nf_inq_varid(ncid3,'med5',varid3)
       if (status.ne.0) write(*,*) "Problem in varid3",status
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read the data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_get_var_real(ncid,varlat,lat)
c       write(*,*) "this is latitude",lat
C
       status = nf_get_var_real(ncid,varlon,lon)
c       write(*,*) "this is longitude",lon
C
       status = nf_get_var_real(ncid,vardep,depth)
       write(*,*) "this is depth",depth
C
c       write(*,*) "this is time",time 
       time(1) = 1998.                   
C
       status = nf_get_var_real(ncid,varid,wflx)              
       write(*,*) "this is chl",(wflx(227,71,1,24))
       status = nf_get_var_int(ncid2,varid2,dmax)              
       write(*,*) "this is dmax",(dmax(227,71,1,nn))
       status = nf_get_var_real(ncid3,varid3,med5)              
       write(*,*) "this is med5",(med5(227,71,1,1))
             
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Re-Locate maximum chl_fill value
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO  ji = 1, jpi
          DO jj = 1, jpj

             IF (dmax(ji,jj,1,nn).NE.fix) then

                date_max(ji,jj) = dmax(ji,jj,1,nn)
                datemaxori(ji,jj) = dmax(ji,jj,1,nn) - jan(nn)

             ELSE

                date_max(ji,jj) = fix
                datemaxori(ji,jj) = fix

             ENDIF

             IF (date_max(ji,jj).NE.fix) then

                maxval(ji,jj) = wflx(ji,jj,1,date_max(ji,jj))            
             ELSE
                maxval(ji,jj) = fix2
             ENDIF

          ENDDO
       ENDDO

       write(*,*) "this is maxval", maxval(227,71)
       write(*,*) "this is week of bloom max val", date_max(227,71)
       write(*,*) "this is week of bloom max val ori", 
     +      datemaxori(227,71)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculate threshold 20% of maxval
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

       DO ji = 1, jpi
          DO jj = 1, jpj

             IF (maxval(ji,jj).GT.0.0) then
                thres(ji,jj) = maxval(ji,jj) - maxval(ji,jj)*0.2
             ELSE
                thres(ji,jj) = fix2
             ENDIF

          ENDDO
       ENDDO
       write(*,*)"this is thres", thres(180,130)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Write med0 matrix of fix values
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO ji = 1, jpi
          DO jj = 1, jpj
             DO it = 1, jit

                med0(ji,jj,it) = fix2
             ENDDO
          ENDDO
       ENDDO

       write(*,*)"this is med0", med0(180,130,2)
  
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculate chl_fill minus threshold
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO ji = 1, jpi
          DO jj = 1, jpj

                IF (date_max(ji,jj).NE.fix) then

                   DO it = tstart(nn), jits(nn)

                      med0(ji,jj,it) = wflx(ji,jj,1,it) - 
     +                     med5(ji,jj,1,1)/1.05*1.2
c     +                     thres(ji,jj)

                   ENDDO

                ENDIF
                   
          ENDDO
       ENDDO
       write(*,*)"this is wflx minus thres",(med0(180,130,it),
     +      it = tstart(nn),jits(nn))

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Locate date_start before date_max
c as the change of sign date from wflxave minus median plus 5%
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO  ji = 1, jpi
          DO jj = 1, jpj

             IF (maxval(ji,jj).GT.0.0.AND.
     +            date_max(ji,jj).NE.fix) then
                
                DO kt = date_max(ji,jj), tstart(nn),-1
                   IF (med0(ji,jj,kt).GT.0..AND.
     +                  med0(ji,jj,kt-1).LT.0.) exit
                ENDDO

                IF (wflx(ji,jj,1,kt-1).EQ.fix2) then
                   date_start(ji,jj) = fix
                   datestartori(ji,jj) = fix
                ELSE
                   
                   date_start(ji,jj) = kt - 1
                   datestartori(ji,jj) = kt - 1 - jan(nn) 
                ENDIF

             ELSE
                date_start(ji,jj) = fix
                datestartori(ji,jj) = fix
             ENDIF
          
          ENDDO
       ENDDO
      write(*,*) "this is date of start of bloom", date_start(227,71)
      write(*,*) "this is datestartori", datestartori(227,71)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Locate date_end after date_max
c as the change of sign date from wflxave minus median plus 5%
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO  ji = 1, jpi
          DO jj = 1, jpj

             IF (maxval(ji,jj).GT.0.0.AND.
     +            date_max(ji,jj).NE.fix) then

                IF(nn.LT.nyr) then
                   DO kt = date_max(ji,jj), jits(nn)
                      IF (med0(ji,jj,kt).GT.0.
     +                     .AND.med0(ji,jj,kt+1).LT.0.) exit
                   ENDDO
                   
                   IF (wflx(ji,jj,1,kt+1).EQ.fix2) then
                      date_end(ji,jj) = fix
                      dateendori(ji,jj) = fix
                   ELSE
                      
                      date_end(ji,jj) = kt + 1
                      dateendori(ji,jj) = kt + 1 - jan(nn) 
                   ENDIF
                   
                  
                ELSEIF(nn.EQ.nyr) then
                   DO kt = date_max(ji,jj), jits(nn)-1
                      IF (med0(ji,jj,kt).GT.0.
     +                     .AND.med0(ji,jj,kt+1).LT.0.) exit
                   ENDDO
                   
                   date_end(ji,jj) = kt + 1
                   dateendori(ji,jj) = kt + 1 - jan(nn) 
                   
                ELSE
                   date_end(ji,jj) = fix
                   dateendori(ji,jj) = fix
                ENDIF

             ELSE
                date_end(ji,jj) = fix
                dateendori(ji,jj) = fix
             ENDIF
          ENDDO
       ENDDO
       write(*,*) "this is date of end of bloom", date_end(227,71) 
       write(*,*) "this is dateendori", dateendori(227,71) 

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Estimate the duration of bloom
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      DO  ji = 1, jpi
         DO jj = 1, jpj
            IF (maxval(ji,jj).GT.0.0) then
               IF (date_end(ji,jj).EQ.fix.OR.date_start(ji,jj)
     +              .EQ.fix) then
                  duration(ji,jj) = fix
               ELSE
                  duration(ji,jj) = date_end(ji,jj)-date_start(ji,jj)+1
               ENDIF
            ELSE
               duration(ji,jj) = fix
            ENDIF
         ENDDO
      ENDDO

      write(*,*) "this is the duration of a bloom", duration(227,71)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c only save the date_start,date_max,date_end 
c if duration is not equal to -1000
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
       DO  ji = 1, jpi
          DO jj = 1, jpj
             IF (duration(ji,jj).LE.0.OR.duration(ji,jj).GT.50) then
                date_start(ji,jj) = fix
                datestartori(ji,jj) = fix
                date_max(ji,jj) = fix
                datemaxori(ji,jj) = fix
                date_end(ji,jj) = fix
                dateendori(ji,jj) = fix
                duration(ji,jj) = fix
             ENDIF
          ENDDO
       ENDDO
       write(*,*) "this is the duration", duration(227,71)
       write(*,*) "this is the date_start",date_start(227,71)
       write(*,*) "this is the datestartori", datestartori(227,71)
       write(*,*) "this is the date_max", date_max(227,71)
       write(*,*) "this is the datemaxori", datemaxori(227,71)
       write(*,*) "this is the date_end", date_end(227,71) 
       write(*,*) "this is the dateendori", dateendori(227,71)    

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  write the grid data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_put_var_real(ncod,vorlon,lon)
      status = nf_put_var_real(ncod,vorlat,lat)
      status = nf_put_var_real(ncod,vordep,depth)
      status = nf_put_var_real(ncod,vortim,time(1))
C                    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  write data to file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_put_var_int(ncod,vorid(1),datestartori)
      status = nf_put_var_int(ncod,vorid(2),datemaxori)
      status = nf_put_var_int(ncod,vorid(3),dateendori)
      status = nf_put_var_int(ncod,vorid(4),duration)      
C                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  close the netCDF file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_close(ncid)
      write(*,*) 'all data read'
      status = nf_close(ncod)
      write(*,*) 'global data written'
C
      ENDDO

      end

