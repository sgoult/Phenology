      program bloomdefGoG_OCCCI_case_CIV
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
      INTEGER jpi,jpj,jpk,jit,nvar,nyr,jitt
      PARAMETER(jpi=481,jpj=361,jpk=1,jit=828,jitt=86,nvar=4)
      PARAMETER(nyr=17)
      INTEGER ji,jj,jk,it,ii,iit,nv,nn
CC
      REAL wflx(jpi,jpj,jpk,jit)
CC
      REAL lat(jpj), lon(jpi), depth(jpk), time(nyr)
      INTEGER fix
      REAL fix2
C
CC      
      INCLUDE 'netcdf.inc'
c
      INTEGER ncid, iret, varid, ncod, vorid(nvar)
      INTEGER varlat, varlon, vardep, vartim
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
      DATA filein/"bdef98.nc","bdef99.nc","bdef00.nc",
     +     "bdef01.nc","bdef02.nc","bdef03.nc","bdef04.nc",
     +     "bdef05.nc","bdef06.nc","bdef07.nc","bdef08.nc",
     +     "bdef09.nc","bdef10.nc","bdef11.nc","bdef12.nc",
     +     "bdef13.nc","bdef14.nc"/
CC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  new variables for sorting
c
      INTEGER last, kt, ptr, first, nit, counter, nber
      INTEGER ngd(jpi,jpj), nbd(jpi,jpj)
      INTEGER datemaxori(jpi,jpj), date_max(jpi,jpj)
      INTEGER date_min(jpi,jpj) , dateminori(jpi,jpj)
      REAL hold
      REAL wflxsort(jpi,jpj,jpk,jit) 
      REAL maxval(jpi,jpj), minval(jpi,jpj)
      DATA nit/46/
c
c For each year
      INTEGER tstart(nyr),jits(nyr),jan(nyr)
      INTEGER tstart2(nyr),jits2(nyr),jan2(nyr)

      DATA jan/46,92,138,184,230,276,322,368,414,460,506,552,598,644,
     +     690,736,782/
c
      DATA tstart/47,93,139,185,231,277,323,369,415,461,507,553,599,
     +     645,691,737,783/
      DATA jits/92,138,184,230,276,322,368,414,460,506,552,598,644,
     +     690,736,782,828/

      DATA jan2/59,105,151,197,243,289,335,381,427,473,519,565,611,
     +     657,703,749,795/
      DATA tstart2/60,106,152,198,244,290,336,382,428,474,520,566,612,
     +     658,704,750,796/
      DATA jits2/90,136,182,228,274,320,366,412,458,504,550,596,642,
     +     688,734,780,826/


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
       status = nf_def_var(ncod,'date_min',nf_int,nvdims,
     +      vardims,vorid(1))
       status = nf_def_var(ncod,'date_max',nf_int,nvdims,
     +      vardims,vorid(2))
       status = nf_def_var(ncod,'minval',nf_float,nvdims,
     +      vardims,vorid(3))
       status = nf_def_var(ncod,'maxval',nf_float,nvdims,
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
       DO nv = 1, 2
       status = nf_put_att_text(ncod,vorid(nv),'units',11,'week')
       status = nf_put_att_int(ncod, vorid(nv),'_FillValue', nf_int,
     +      1,fix)
       ENDDO
       DO nv = 3, 4
       status = nf_put_att_text(ncod,vorid(nv),'units',11,'mgChl/m^3')
       status = nf_put_att_real(ncod, vorid(nv),'_FillValue', nf_float,
     +      1,fix2)
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
     +      "OCCCI_v3-8DAY-97_14_9km_pheno_fav8_sbx3weeks.nc"
     +      ,nf_nowrite,ncid)
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Locate maximum chl_fill value
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO  ji = 1, jpi
          DO jj = 1, jpj

             ptr = tstart2(nn)
             first = tstart2(nn) + 1
             DO kt = first, jits2(nn)

                IF (wflx(ji,jj,1,kt).gt.
     +               wflx(ji,jj,1,ptr)) then
                   ptr=kt
                ENDIF
             
             ENDDO
             
             maxval(ji,jj) = wflx(ji,jj,1,ptr)

             IF (maxval(ji,jj).LE.0.0) then
                date_max(ji,jj) = fix
                datemaxori(ji,jj) = fix
             ELSE
                date_max(ji,jj) = ptr
                datemaxori(ji,jj) = ptr - jan(nn)
             ENDIF

          ENDDO
       ENDDO

       write(*,*) "this is maxval", maxval(227,71)
       write(*,*) "this is week of bloom max val", date_max(227,71)
       write(*,*) "this is week of bloom max val ori", 
     +      datemaxori(227,71)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  save wflxsort
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO ji = 1,jpi
          DO jj = 1, jpj
             DO jk = 1, jpk
                DO it = tstart(nn), jits(nn)
                   wflxsort(ji,jj,jk,it) = wflx(ji,jj,jk,it)
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       write(*,*) "this wflxsort", (wflxsort(227,71,1,it),
     +      it = tstart(nn), jits(nn))
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Sort the chl data (for one year)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO  ji = 1, jpi
         DO jj = 1, jpj

               last = jits(nn)
               DO it = tstart(nn), jits(nn) - 1
                  ptr = it
                  first = it + 1
          
                  DO kt = first, last
             
                     IF (ABS(wflxsort(ji,jj,1,kt)).lt.
     +                    ABS(wflxsort(ji,jj,1,ptr))) then
                        ptr=kt
                     ENDIF
             
                  ENDDO
         
                  hold = wflxsort(ji,jj,1,it)
                  wflxsort(ji,jj,1,it) = wflxsort(ji,jj,1,ptr)
                  wflxsort(ji,jj,1,ptr) = hold
 
               ENDDO
         ENDDO
       ENDDO

   
       DO  ji = 1, jpi
          DO jj = 1, jpj
             counter =  0
             DO it = tstart(nn),jits(nn)
                IF (ABS(wflxsort(ji,jj,1,it)).lt.0.5E20) then
                   counter=counter+1
                ENDIF
             ENDDO  
             nbd(ji,jj) = nit - counter
             ngd(ji,jj) = counter
          ENDDO
       ENDDO

       write(*,*) "this is nbd",nbd(227,71)
       write(*,*) "this is ngd",ngd(227,71)
       write(*,*) "this is sorted chl",
     +    (wflxsort(227,71,1,it),it = tstart(nn), jits(nn))
       write(*,*) "this is sorted good chl data",
     +    (wflxsort(227,71,1,it),it = tstart(nn), tstart(nn) + 
     +      ngd(227,71)-1)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Find minimum chl_fill value
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO  ji = 1, jpi
          DO jj = 1, jpj

             ptr = tstart(nn)
             first = tstart(nn) + 1
             DO kt = first, ngd(ji,jj)
                
                IF (wflxsort(ji,jj,1,kt).lt.
     +               wflxsort(ji,jj,1,ptr)) then
                   ptr=kt
                ENDIF
                
             ENDDO
             
             minval(ji,jj) = wflxsort(ji,jj,1,ptr)
             
             
          ENDDO
       ENDDO

       write(*,*) "this is minval", minval(180,130)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Locate minimum chl_fill value
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO ji = 1, jpi
          DO jj = 1, jpj

             ptr = tstart(nn)
             first = tstart(nn)+1
             DO kt = first, jits(nn)

                IF (wflx(ji,jj,1,kt).EQ.minval(ji,jj)) then
                   ptr=kt
                ENDIF
             
             ENDDO
             
             IF (minval(ji,jj).GE.0.0) then
                date_min(ji,jj) = ptr
                dateminori(ji,jj) = ptr - jan(nn)
             ELSE
                date_min(ji,jj) = fix 
                dateminori(ji,jj) = fix
             ENDIF
             
          ENDDO
       ENDDO

       write(*,*) "this is minval", minval(180,130)
       write(*,*) "this is week of bloom min val", date_min(180,130)
c       write(*,*) "thisis date_min", (date_min(ji,130), ji = 1,jpi) 
c       write(*,*) "thisis date_min", (date_min(180,jj), jj = 1,jpj)
       write(*,*) "this is week of bloom min val ori", 
     +      dateminori(180,130)


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
      status = nf_put_var_int(ncod,vorid(1),date_min)
      status = nf_put_var_int(ncod,vorid(2),date_max)
      status = nf_put_var_real(ncod,vorid(3),minval)
      status = nf_put_var_real(ncod,vorid(4),maxval)      
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

