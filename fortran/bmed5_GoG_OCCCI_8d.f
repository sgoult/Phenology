      program bloom_median_GoG_OCCCI_8d
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
      PARAMETER(jpi=481,jpj=361,jpk=1,jit=798,nvar=4)
      PARAMETER(nyr=18)
      INTEGER ji,jj,jk,it,ii,iit,nv,nn
CC
      REAL wflx(jpi,jpj,jpk,jit)
CC
      REAL lat(jpj), lon(jpi), depth(jpk), time(jit)
      INTEGER fix
      REAL fix2
C
CC      
      INCLUDE 'netcdf.inc'
c
      INTEGER ncid, iret, varid, ncod, vorid
      INTEGER ncid2, varid2(nvar)
      INTEGER varlat, varlon, vardep, vartim
      INTEGER varlat2, varlon2, vardep2, vartim2
      INTEGER vorlat, vorlon, vordep, vortim
      INTEGER nvdims
      PARAMETER(nvdims = 4)
      INTEGER vardims(nvdims),start(nvdims),count(nvdims),stride(nvdims)
      INTEGER status

CC
      DATA start/1,1,1,1/
      DATA count/jpi,jpj,jpk,jit/
      DATA stride/1,1,1,1/
      DATA fix,fix2/-1000,-1.e+34/
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  new variables for sorting
c
      INTEGER last, kt, ptr, first, nit, counter, nber
      INTEGER ngd(jpi,jpj), nbd(jpi,jpj)
      REAL median(jpi,jpj), med5(jpi,jpj)
      REAL wflxsort(jpi,jpj,jpk,jit) , hold
      DATA nit/46/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  open the new netCDF file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
       status = nf_create("bmed5_GoG_OCCCI_8d.nc", nf_clobber,ncod)
c       
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
       status = nf_def_var(ncod,'med5',nf_float,nvdims,vardims,vorid)
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
       status = nf_put_att_text(ncod,vorid,'units',11,'mgChl/m3')
       status = nf_put_att_real(ncod, vorid,'_FillValue', nf_float,
     +      1,fix2)
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
     +"CCI_ALL-v3.0-8DAY-Sep97-Dec14_9kmba.nc",nf_nowrite,ncid)
      write(*,*) "File open", status
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  inquire attributes for the grid and variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_inq_varid(ncid,'LONGITUDE6181_6661',varlon)
       if (status.ne.0) write(*,*) "Problem in varlon",status
C
       status = nf_inq_varid(ncid,'LATITUDE900_1260',varlat)
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
c  save wflxsort
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO ji = 1,jpi
          DO jj = 1, jpj
             DO jk = 1, jpk
                DO it = 1, jit
                   wflxsort(ji,jj,jk,it) = wflx(ji,jj,jk,it)
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       write(*,*) "this wflxsort", (wflxsort(227,71,1,it),
     +      it = 1, jit)            
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Sort the chl data over the 18 years
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO  ji = 1, jpi
         DO jj = 1, jpj

               last = jit
               DO it = 1, jit - 1
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
             DO it = 1,jit
                IF (ABS(wflxsort(ji,jj,1,it)).lt.0.5E20) then
                   counter=counter+1
                ENDIF
             ENDDO  
             nbd(ji,jj) = jit - counter
             ngd(ji,jj) = counter
          ENDDO
       ENDDO

       write(*,*) "this is nbd",nbd(227,71)
       write(*,*) "this is ngd",ngd(227,71)
       write(*,*) "this is sorted chl",
     +    (wflxsort(227,71,1,it),it = 1, jit)
       write(*,*) "this is sorted good chl data",
     +    (wflxsort(227,71,1,it),it = 1, 1+ngd(227,71)-1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Find median value over the 18 years
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO ji = 1, jpi
          DO jj = 1, jpj
             IF (MOD(ngd(ji,jj),2).ne.0) then
                median(ji,jj) = wflxsort(ji,jj,1,
     +               ngd(ji,jj)/2+1)
             ELSE 
                median(ji,jj) = (wflxsort(ji,jj,1,
     +               ngd(ji,jj)/2)+
     +            wflxsort(ji,jj,1,1-1+ngd(ji,jj)/2+1))/2.0
             ENDIF
          ENDDO
       ENDDO
             
       write(*,*)"this is median value", median(227,71)  
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Median value plus 5% 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       DO ji = 1, jpi
          DO jj = 1, jpj
             IF (median(ji,jj).GE.0.0) then
                med5(ji,jj) = median(ji,jj) + median(ji,jj)*0.05
             ELSE
                med5(ji,jj) = fix2
             ENDIF
          ENDDO
       ENDDO
       write(*,*)"this is median plus 5%", med5(227,71) 
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
      status = nf_put_var_real(ncod,vorid,med5) 
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

      end

