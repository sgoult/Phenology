      program save_chl_file_jan97_dec14
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
      INTEGER jpi,jpj,jpk,jit,jit2
      PARAMETER(jpi=481,jpj=361,jpk=1,jit=844,jit2=828)
      INTEGER ji,jj,jk,it,jt
CC
      CHARACTER*52 namein
CC
      REAL wflxchl(jpi,jpj,jpk,jit)
CC
      REAL lat(jpj), lon(jpi), depth(jpk), time(jit2)
      REAL fix
CC
      INCLUDE 'netcdf.inc'
c
      INTEGER ncid, iret, varid, ncod, vorid
      INTEGER varlat, varlon, vardep, vartim
      INTEGER vorlat, vorlon, vordep, vortim
      INTEGER nvdims
      PARAMETER(nvdims = 4)
      INTEGER vardims(nvdims),start(nvdims),count(nvdims),stride(nvdims)
      INTEGER vordims(nvdims)
      INTEGER status
c
      DATA start/1,1,1,1/
      DATA count/jpi,jpj,jpk,jit2/
      DATA stride/1,1,1,1/
      DATA fix/-1.e+34/

      REAL wflxchlt(jpi,jpj,jpk,jit2)

C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  open the new netCDF file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      status = nf_create(
     +     'CCIv3-8DAY-Jan97-Dec14_9km_fillphenob.nc', 
     +     nf_clobber,ncod)
C       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define attributes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_def_dim(ncod,'LONGITUDE',jpi,vardims(1))
       status = nf_def_dim(ncod,'LATITUDE',jpj,vardims(2))
       status = nf_def_dim(ncod,'DEPTH',jpk,vardims(3))
       status = nf_def_dim(ncod,'TIME',jit2,vardims(4))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  
       status=nf_def_var(ncod,'chl_fill',nf_float,nvdims,vardims,vorid)
       status=nf_def_var(ncod,'LONGITUDE',nf_float,1,vardims(1),vorlon)
       status=nf_def_var(ncod,'LATITUDE',nf_float,1,vardims(2),vorlat)
       status=nf_def_var(ncod,'DEPTH',nf_float,1,vardims(3),vordep)
       status=nf_def_var(ncod,'TIME',nf_float,1,vardims(4),vortim)
C   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  create attributes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
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
       status = nf_put_att_text(ncod, vorid, 'units', 11, 'mg/m3')
       status = nf_put_att_real(ncod, vorid, '_FillValue', nf_float,1,
     +                                                           fix)
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
       status = nf_open(
     +      "CCI_ALL-v3.0-8DAY-Sep97-Dec15_9kmba_fillphenob.nc", 
     +      nf_nowrite,ncid)
      write(*,*) "File open chl", nf_strerror(status)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  inquire attributes for the grid and variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_inq_varid(ncid,'LONGITUDE6181_6661',varlon)
       if (status.ne.0) write(*,*) "Problem in varlon",
     +      nf_strerror(status)
C
       status = nf_inq_varid(ncid,'LATITUDE900_1260',varlat)
       if (status.ne.0) write(*,*) "Problem in varlat",status
C
c       status = nf_inq_varid(ncid,'DEPTH',vardep)
c       if (status.ne.0) write(*,*) "Problem in vardep",status
C
c       status = nf_inq_varid(ncid,'TIME',vartim)
c       if (status.ne.0) write(*,*) "Problem in vartim",status
C
       status = nf_inq_varid(ncid,'CHL_FILL',varid)
       if (status.ne.0) write(*,*) "Problem in varid",status
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read the data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_get_var_real(ncid,varlat,lat)
       write(*,*) 'get lat', nf_strerror(status)
       write(*,*) "this is latitude",lat(263)
c       write(*,*) "this is latitude0",lat(881)
C
       status = nf_get_var_real(ncid,varlon,lon)
       write(*,*) 'get lon', nf_strerror(status)
       write(*,*) "this is longitude",lon(12)
C
c       status = nf_get_var_real(ncid,vardep,depth)
c       write(*,*) "this is depth",depth
C
c       status = nf_get_var_real(ncid,vartim,time)
c       write(*,*) "this is time",time                    
C
       status = nf_get_var_real(ncid,varid,wflxchl)              
       write(*,*) "this is wflxchl",wflxchl(220,300,1,184)    
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  save data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
       DO ji = 1, jpi
          DO jj = 1, jpj
             DO it = 1, 30
                wflxchlt(ji,jj,1,it) = fix
             ENDDO

             DO it = 31, jit2
                wflxchlt(ji,jj,1,it) = wflxchl(ji,jj,1,it-30)
             ENDDO

          ENDDO
       ENDDO
       
       write(*,*) "this is wflxchlt", wflxchlt(22,71,1,84) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  save time
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
       DO it = 1, jit2
          time(it) = real(it)
       ENDDO

       depth(1) = 0.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  write the grid data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_put_var_real(ncod,vorlon,lon)
      status = nf_put_var_real(ncod,vorlat,lat)
      status = nf_put_var_real(ncod,vordep,depth)
      status = nf_put_var_real(ncod,vortim,time)
C                    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  write data to file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_put_var_real(ncod,vorid,wflxchlt)
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

