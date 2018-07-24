      program mergebdefGoG_OCCCI_c2_CIV
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
      INTEGER jpi,jpj,jpk,jit,nvar
      PARAMETER(jpi=481,jpj=361,jpk=1,jit=17,nvar=4)
      INTEGER ji,jj,jk,it,nv,nn
CC
      REAL wflx(jpi,jpj)
CC
      REAL lat(jpj), lon(jpi), depth(jpk), time(jit)
      INTEGER fix
      REAL fix2
CC
      INCLUDE 'netcdf.inc'
c
      INTEGER ncid, iret, varid(nvar), ncod,vorid(nvar)
      INTEGER ncid1,ncid2,ncid3,ncid4,ncid5,ncid6,ncid7,ncid8 
      INTEGER ncid9,ncid10,ncid11,ncid12,ncid13,ncid14,ncid15,ncid16 
      INTEGER varid1(nvar),varid2(nvar),varid3(nvar),varid4(nvar)
      INTEGER varid5(nvar),varid6(nvar),varid7(nvar),varid8(nvar)
      INTEGER varid9(nvar),varid10(nvar),varid11(nvar),varid12(nvar)
      INTEGER varid13(nvar),varid14(nvar),varid15(nvar),varid16(nvar)  
      INTEGER ncid17,varid17(nvar)
      INTEGER varlat, varlon, vardep, vartim, vartim2
      INTEGER vorlat, vorlon, vordep, vortim
      INTEGER nvdims
      PARAMETER(nvdims = 4)
      INTEGER vardims(nvdims),start(nvdims),count(nvdims),stride(nvdims)
      INTEGER status
      INTEGER date_start(jpi,jpj), date_max(jpi,jpj), date_end(jpi,jpj)
      INTEGER duration(jpi,jpj)
      CHARACTER*20 varname(nvar)
c
      DATA start/1,1,1,1/
      DATA count/jpi,jpj,jpk,1/
      DATA stride/1,1,1,1/
      DATA fix,fix2/-1000,-1.e+34/

      DATA varname/'date_start','date_max','date_end','duration'/

CC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  open the netCDF file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_create("mergedbdefGoG_OCCCI_c2_CIV_v4.nc",
     +     nf_clobber,ncod)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define attributes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_def_dim(ncod,'LONGITUDE',jpi,vardims(1))
       status = nf_def_dim(ncod,'LATITUDE',jpj,vardims(2))
       status = nf_def_dim(ncod,'DEPTH',jpk,vardims(3))
       status = nf_def_dim(ncod,'TIME',jit,vardims(4)) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  
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
       DO nv = 1, nvar
       status = nf_put_att_text(ncod,vorid(nv),'units',11,'weeks')
       status =nf_put_att_int(ncod, vorid(nv),'_FillValue', nf_int,1,
     +                                                           fix)
       ENDDO  
c
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  end of define mode
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_enddef(ncod)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  open the netCDF input file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C NOTE: the variable status should always return a 0 value. When this
C is not 0, then the operation failed (check the name of the file and
C variables, watch for capital letters)
Cbdateendnorth_00_01.nc

       status = nf_open("bddef98.nc",nf_nowrite,ncid)
       status = nf_open("bddef99.nc",nf_nowrite,ncid2)
       status = nf_open("bddef00.nc",nf_nowrite,ncid3)
       status = nf_open("bddef01.nc",nf_nowrite,ncid4)
       status = nf_open("bddef02.nc",nf_nowrite,ncid5)
       status = nf_open("bddef03.nc",nf_nowrite,ncid6)
       status = nf_open("bddef04.nc",nf_nowrite,ncid7)
       status = nf_open("bddef05.nc",nf_nowrite,ncid8)
       status = nf_open("bddef06.nc",nf_nowrite,ncid9)
       status = nf_open("bddef07.nc",nf_nowrite,ncid10)
       status = nf_open("bddef08.nc",nf_nowrite,ncid11)
       status = nf_open("bddef09.nc",nf_nowrite,ncid12)
       status = nf_open("bddef10.nc",nf_nowrite,ncid13)
       status = nf_open("bddef11.nc",nf_nowrite,ncid14)
       status = nf_open("bddef12.nc",nf_nowrite,ncid15)
       status = nf_open("bddef13.nc",nf_nowrite,ncid16)  
       status = nf_open("bddef14.nc",nf_nowrite,ncid17)  

       write(*,*) "File open", status
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
C      status = nf_inq_varid(ncid,'TIME',vartim)
C      status = nf_inq_varid(ncid2,'TIME',vartim2)
C      if (status.ne.0) write(*,*) "Problem in vartim",status
C
       DO nv = 1, nvar
          status = nf_inq_varid(ncid,varname(nv),varid(nv))
          status = nf_inq_varid(ncid2,varname(nv),varid2(nv))
          status = nf_inq_varid(ncid3,varname(nv),varid3(nv))
          status = nf_inq_varid(ncid4,varname(nv),varid4(nv))
          status = nf_inq_varid(ncid5,varname(nv),varid5(nv))
          status = nf_inq_varid(ncid6,varname(nv),varid6(nv))
          status = nf_inq_varid(ncid7,varname(nv),varid7(nv))
          status = nf_inq_varid(ncid8,varname(nv),varid8(nv))
          status = nf_inq_varid(ncid9,varname(nv),varid9(nv))
          status = nf_inq_varid(ncid10,varname(nv),varid10(nv))
          status = nf_inq_varid(ncid11,varname(nv),varid11(nv))
          status = nf_inq_varid(ncid12,varname(nv),varid12(nv))
          status = nf_inq_varid(ncid13,varname(nv),varid13(nv))
          status = nf_inq_varid(ncid14,varname(nv),varid14(nv))
          status = nf_inq_varid(ncid15,varname(nv),varid15(nv))
          status = nf_inq_varid(ncid16,varname(nv),varid16(nv))       
          status = nf_inq_varid(ncid17,varname(nv),varid17(nv))  
       ENDDO

       if (status.ne.0) write(*,*) "Problem in varid",status
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read the grid data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       status = nf_get_var_real(ncid,varlat,lat)
       write(*,*) "this is latitude",lat
C
       status = nf_get_var_real(ncid,varlon,lon)
       write(*,*) "this is longitude",lon
C
       status = nf_get_var_real(ncid,vardep,depth)
       write(*,*) "this is depth",depth
                
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  write the grid data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_put_var_real(ncod,vorlon,lon)
      status = nf_put_var_real(ncod,vorlat,lat)
      status = nf_put_var_real(ncod,vordep,depth)
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read and write the data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C do 1998
C
      status = nf_get_vara_int(ncid,varid(1),start,count,date_start)
      status = nf_get_vara_int(ncid,varid(2),start,count,date_max)
      status = nf_get_vara_int(ncid,varid(3),start,count,date_end)
      status = nf_get_vara_int(ncid,varid(4),start,count,duration)
      time(1) = 1998.0
      write(*,*) "test1, date_start ",date_start(60,60)
      write(*,*) "test1, date_max ",date_max(60,60)
      write(*,*) "test1, date_end ",date_end(60,60)
      write(*,*) "test1, duration ",duration(60,60)
  
      start(4) = 1
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 1999
C
      status = nf_get_vara_int(ncid2,varid2(1),start,count,date_start)
      status = nf_get_vara_int(ncid2,varid2(2),start,count,date_max)
      status = nf_get_vara_int(ncid2,varid2(3),start,count,date_end)
      status = nf_get_vara_int(ncid2,varid2(4),start,count,duration)
      time(2) = 1999.0
      write(*,*) "test2, date_start ",date_start(60,60)
      write(*,*) "test2, date_max ",date_max(60,60)
      write(*,*) "test2, date_end ",date_end(60,60)
      write(*,*) "test2, duration ",duration(60,60)

      start(4) = 2
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2000
C
      start(4) = 1
      status = nf_get_vara_int(ncid3,varid3(1),start,count,date_start)
      status = nf_get_vara_int(ncid3,varid3(2),start,count,date_max)
      status = nf_get_vara_int(ncid3,varid3(3),start,count,date_end)
      status = nf_get_vara_int(ncid3,varid3(4),start,count,duration) 
      time(3) = 2000.0
      write(*,*) "test3, date_start ",date_start(60,60)
      write(*,*) "test3, date_max ",date_max(60,60)
      write(*,*) "test3, date_end ",date_end(60,60)
      write(*,*) "test3, duration ",duration(60,60)

      start(4) = 3
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2001
C
      start(4) = 1
      status = nf_get_vara_int(ncid4,varid4(1),start,count,date_start)
      status = nf_get_vara_int(ncid4,varid4(2),start,count,date_max)
      status = nf_get_vara_int(ncid4,varid4(3),start,count,date_end)
      status = nf_get_vara_int(ncid4,varid4(4),start,count,duration) 
      time(4) = 2001.0
      write(*,*) "test4, date_start ",date_start(60,60)
      write(*,*) "test4, date_max ",date_max(60,60)
      write(*,*) "test4, date_end ",date_end(60,60)
      write(*,*) "test4, duration ",duration(60,60)  

      start(4) = 4
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
c
C do 2002
C
      start(4) = 1
      status = nf_get_vara_int(ncid5,varid5(1),start,count,date_start)
      status = nf_get_vara_int(ncid5,varid5(2),start,count,date_max)
      status = nf_get_vara_int(ncid5,varid5(3),start,count,date_end)
      status = nf_get_vara_int(ncid5,varid5(4),start,count,duration)   
        
      time(5) = 2002.0
      write(*,*) "test5, date_start ",date_start(60,60)
      write(*,*) "test5, date_max ",date_max(60,60)
      write(*,*) "test5, date_end ",date_end(60,60)
      write(*,*) "test5, duration ",duration(60,60)  

      start(4) = 5
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
   
C do 2003
C
      start(4) = 1
      status = nf_get_vara_int(ncid6,varid6(1),start,count,date_start)
      status = nf_get_vara_int(ncid6,varid6(2),start,count,date_max)
      status = nf_get_vara_int(ncid6,varid6(3),start,count,date_end)
      status = nf_get_vara_int(ncid6,varid6(4),start,count,duration)
   
      time(6) = 2003.0
      write(*,*) "test6, date_start ",date_start(60,60)
      write(*,*) "test6, date_max ",date_max(60,60)
      write(*,*) "test6, date_end ",date_end(60,60)
      write(*,*) "test6, duration ",duration(60,60)  

      start(4) = 6
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2004
C
      start(4) = 1
      status = nf_get_vara_int(ncid7,varid7(1),start,count,date_start)
      status = nf_get_vara_int(ncid7,varid7(2),start,count,date_max)
      status = nf_get_vara_int(ncid7,varid7(3),start,count,date_end)
      status = nf_get_vara_int(ncid7,varid7(4),start,count,duration)
           
      time(7) = 2004.0
      write(*,*) "test7, date_start ",date_start(60,60)
      write(*,*) "test7, date_max ",date_max(60,60)
      write(*,*) "test7, date_end ",date_end(60,60)
      write(*,*) "test7, duration ",duration(60,60) 

      start(4) = 7
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2005
C
      start(4) = 1
      status = nf_get_vara_int(ncid8,varid8(1),start,count,date_start)
      status = nf_get_vara_int(ncid8,varid8(2),start,count,date_max)
      status = nf_get_vara_int(ncid8,varid8(3),start,count,date_end)
      status = nf_get_vara_int(ncid8,varid8(4),start,count,duration)
           
      time(8) = 2005.0
      write(*,*) "test8, date_start ",date_start(60,60)
      write(*,*) "test8, date_max ",date_max(60,60)
      write(*,*) "test8, date_end ",date_end(60,60)
      write(*,*) "test8, duration ",duration(60,60) 

      start(4) = 8
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2006
C
      start(4) = 1
      status = nf_get_vara_int(ncid9,varid9(1),start,count,date_start)
      status = nf_get_vara_int(ncid9,varid9(2),start,count,date_max)
      status = nf_get_vara_int(ncid9,varid9(3),start,count,date_end)
      status = nf_get_vara_int(ncid9,varid9(4),start,count,duration)
           
      time(9) = 2006.0
      write(*,*) "test9, date_start ",date_start(60,60)
      write(*,*) "test9, date_max ",date_max(60,60)
      write(*,*) "test9, date_end ",date_end(60,60)
      write(*,*) "test9, duration ",duration(60,60) 

      start(4) = 9
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2007
C
      start(4) = 1      
      status = nf_get_vara_int(ncid10,varid10(1),start,count,date_start)
      status = nf_get_vara_int(ncid10,varid10(2),start,count,date_max)
      status = nf_get_vara_int(ncid10,varid10(3),start,count,date_end)
      status = nf_get_vara_int(ncid10,varid10(4),start,count,duration)
      time(10) = 2007.0
      write(*,*) "test2, date_start ",date_start(60,60)
      write(*,*) "test2, date_max ",date_max(60,60)
      write(*,*) "test2, date_end ",date_end(60,60)
      write(*,*) "test2, duration ",duration(60,60)

      start(4) = 10
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2008
C
      start(4) = 1
      status = nf_get_vara_int(ncid11,varid11(1),start,count,date_start)
      status = nf_get_vara_int(ncid11,varid11(2),start,count,date_max)
      status = nf_get_vara_int(ncid11,varid11(3),start,count,date_end)
      status = nf_get_vara_int(ncid11,varid11(4),start,count,duration) 
      time(11) = 2008.0
      write(*,*) "test3, date_start ",date_start(60,60)
      write(*,*) "test3, date_max ",date_max(60,60)
      write(*,*) "test3, date_end ",date_end(60,60)
      write(*,*) "test3, duration ",duration(60,60)

      start(4) = 11
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2009
C
      start(4) = 1
      status = nf_get_vara_int(ncid12,varid12(1),start,count,date_start)
      status = nf_get_vara_int(ncid12,varid12(2),start,count,date_max)
      status = nf_get_vara_int(ncid12,varid12(3),start,count,date_end)
      status = nf_get_vara_int(ncid12,varid12(4),start,count,duration) 
      time(12) = 2009.0
      write(*,*) "test4, date_start ",date_start(60,60)
      write(*,*) "test4, date_max ",date_max(60,60)
      write(*,*) "test4, date_end ",date_end(60,60)
      write(*,*) "test4, duration ",duration(60,60)  

      start(4) = 12
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
c
C do 2010
C
      start(4) = 1
      status = nf_get_vara_int(ncid13,varid13(1),start,count,date_start)
      status = nf_get_vara_int(ncid13,varid13(2),start,count,date_max)
      status = nf_get_vara_int(ncid13,varid13(3),start,count,date_end)
      status = nf_get_vara_int(ncid13,varid13(4),start,count,duration)   
        
      time(13) = 2010.0
      write(*,*) "test5, date_start ",date_start(60,60)
      write(*,*) "test5, date_max ",date_max(60,60)
      write(*,*) "test5, date_end ",date_end(60,60)
      write(*,*) "test5, duration ",duration(60,60)  

      start(4) = 13
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
   
C do 2011
C
      start(4) = 1
      status = nf_get_vara_int(ncid14,varid14(1),start,count,date_start)
      status = nf_get_vara_int(ncid14,varid14(2),start,count,date_max)
      status = nf_get_vara_int(ncid14,varid14(3),start,count,date_end)
      status = nf_get_vara_int(ncid14,varid14(4),start,count,duration)
   
      time(14) = 2011.0
      write(*,*) "test6, date_start ",date_start(60,60)
      write(*,*) "test6, date_max ",date_max(60,60)
      write(*,*) "test6, date_end ",date_end(60,60)
      write(*,*) "test6, duration ",duration(60,60)  

      start(4) = 14
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2012
C
      start(4) = 1
      status = nf_get_vara_int(ncid15,varid15(1),start,count,date_start)
      status = nf_get_vara_int(ncid15,varid15(2),start,count,date_max)
      status = nf_get_vara_int(ncid15,varid15(3),start,count,date_end)
      status = nf_get_vara_int(ncid15,varid15(4),start,count,duration)
           
      time(15) = 2012.0
      write(*,*) "test7, date_start ",date_start(60,60)
      write(*,*) "test7, date_max ",date_max(60,60)
      write(*,*) "test7, date_end ",date_end(60,60)
      write(*,*) "test7, duration ",duration(60,60) 

      start(4) = 15
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2013
C
      start(4) = 1
      status = nf_get_vara_int(ncid16,varid16(1),start,count,date_start)
      status = nf_get_vara_int(ncid16,varid16(2),start,count,date_max)
      status = nf_get_vara_int(ncid16,varid16(3),start,count,date_end)
      status = nf_get_vara_int(ncid16,varid16(4),start,count,duration)
           
      time(16) = 2013.0
      write(*,*) "test8, date_start ",date_start(60,60)
      write(*,*) "test8, date_max ",date_max(60,60)
      write(*,*) "test8, date_end ",date_end(60,60)
      write(*,*) "test8, duration ",duration(60,60) 

      start(4) = 16
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
c
C do 2014
C
      start(4) = 1
      status = nf_get_vara_int(ncid17,varid17(1),start,count,date_start)
      status = nf_get_vara_int(ncid17,varid17(2),start,count,date_max)
      status = nf_get_vara_int(ncid17,varid17(3),start,count,date_end)
      status = nf_get_vara_int(ncid17,varid17(4),start,count,duration)
           
      time(17) = 2014.0
      write(*,*) "test17, date_start ",date_start(60,60)
      write(*,*) "test17, date_max ",date_max(60,60)
      write(*,*) "test17, date_end ",date_end(60,60)
      write(*,*) "test17, duration ",duration(60,60) 

      start(4) = 17
      status = nf_put_vara_int(ncod,vorid(1),start,count,date_start)
      status = nf_put_vara_int(ncod,vorid(2),start,count,date_max)
      status = nf_put_vara_int(ncod,vorid(3),start,count,date_end)
      status = nf_put_vara_int(ncod,vorid(4),start,count,duration)
C
C write time axis
C
      status = nf_put_var_real(ncod,vortim,time)
C
      write(*,*) "data written for ", time
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  close the netCDF file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      status = nf_close(ncid)
      status = nf_close(ncid2)
      status = nf_close(ncid3)
      status = nf_close(ncid4)
      status = nf_close(ncid5)
      status = nf_close(ncid6)
      status = nf_close(ncid7)
      status = nf_close(ncid8)
      status = nf_close(ncid9)
      status = nf_close(ncid10)
      status = nf_close(ncid11)
      status = nf_close(ncid12)
      status = nf_close(ncid13)
      status = nf_close(ncid14)
      status = nf_close(ncid15)
      status = nf_close(ncid16)
      status = nf_close(ncid17)
      status = nf_close(ncod)
      write(*,*) 'all data read'
C
      END

