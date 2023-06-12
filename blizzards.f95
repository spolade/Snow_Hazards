program blizzards
implicit none
 
! The program calculates the total amount of blizzard days during the period
! xxxx-xxxx in each grid cell, and on how many years blizzard days occur
!
! Mean wind speed is sometimes used instead of maximum wind gusts
! The precipitation levels in the MPI model need to be multiplied by four!
 
integer,parameter :: ii=221, jj=103, day_max=360! day_max=366
integer,parameter :: yr=30
integer :: irec, day, i, j, ipx, jpx, yrbegin, yrend, y, yrlength, yq
real,dimension(ii,jj,day_max,yr) :: t2m, prec, gust
real,dimension(ii,jj,yr) :: cases
real,dimension(ii,jj) :: cases_tot, cases_annual
real :: lim_t2m=0, lim_prec=10, lim_gust=17! lim_gust=7
 
 
!an example grid cell
ipx=138 !lon=23.75E
jpx=79 !lat=60.75N
 
yrbegin=1971
yrend=yrbegin+yr-1
 
! Read in the daily mean temperatures
open(11,file='/RAIN/temporal/blizzards/&
tas_EUR-44reg_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22E_v1_day_1971-2000.gr',&
     form='unformatted',status='old',access='direct',recl=4*ii*jj)
 
irec=1
 
do y=yrbegin,yrend
   if(mod(y,4)==0.and.day_max==366)then
      if(y/=2100)then
         yrlength=366
      end if
   else
      yrlength=365
   end if
   if(day_max==360)then
      yrlength=day_max
      if(y==yrend)then
         !yrlength=330
      end if
   end if
   yq=y-yrbegin+1
   do day=1,yrlength
      read(11,rec=irec)((t2m(i,j,day,yq),i=1,ii),j=1,jj)
      irec=irec+1
   end do
end do
 
close (11)
 
do i=1,ii
   do j=1,jj
      do y=yrbegin,yrend
         if(mod(y,4)==0.and.day_max==366)then
            if(v/=2100)then
               yrlength=366
            end if
         else
            yrlength=365
         end if
         if(day_max==360)then
            yrlength=day_max
            if(y==yrend)then
               !yrlength=330
            end if
         end if
         yq=y-yrbegin+1
         do day=1,yrlength
            t2m(i,j,day,yq)=t2m(i,j,day,yq)-273.15
         end do
      end do
   end do
end do
 
! Read in the daily precipitation levels
open(11,file='/lustre/tmp/lehtonei/RAIN/temporal/blizzards/&
pr_EUR-44reg_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22E_v1_day_1971-2000.gr',&
     form='unformatted',status='old',access='direct',recl=4*ii*jj)
 
irec=1
 
do y=yrbegin,yrend
   if(mod(y,4)==0.and.day_max==366)then
      if(v/=2100)then
         yrlength=366
      end if
   else
      yrlength=365
   end if
   if(day_max==360)then
      yrlength=day_max
      if(y==yrend)then
         !yrlength=330
      end if
   end if
   yq=y-yrbegin+1
   do day=1,yrlength
      read(11,rec=irec)((prec(i,j,day,yq),i=1,ii),j=1,jj)
      irec=irec+1
   end do
end do
 
close (11)
 
do i=1,ii
   do j=1,jj
      do y=yrbegin,yrend
         if(mod(y,4)==0.and.day_max==366)then
            if(v/=2100)then
               yrlength=366
            end if
         else
            yrlength=365
         end if
         if(day_max==360)then
            yrlength=day_max
            if(y==yrend)then
               !yrlength=330
            end if
         end if
         yq=y-yrbegin+1
         do day=1,yrlength
            prec(i,j,day,yq)=prec(i,j,day,yq)*86400
!            prec(i,j,day,yq)=prec(i,j,day,yq)*4 !For the MPI model
         end do
      end do
   end do
end do
 
! Read in the winds (gusts)
open(11,file='/lustre/tmp/lehtonei/RAIN/temporal/blizzards/&
wsgsmax_EUR-44reg_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22E_v1_day_1971-2000.gr',&
     form='unformatted',status='old',access='direct',recl=4*ii*jj)
 
irec=1
 
do y=yrbegin,yrend
   if(mod(y,4)==0.and.day_max==366)then
      if(v/=2100)then
         yrlength=366
      end if
   else
      yrlength=365
   end if
   if(day_max==360)then
      yrlength=day_max
      if(y==yrend)then
         !yrlength=330
      end if
   end if
   yq=y-yrbegin+1
   do day=1,yrlength
      read(11,rec=irec)((gust(i,j,day,yq),i=1,ii),j=1,jj)
      irec=irec+1
   end do
end do
 
close (11)
 
! Calculate the numer of blizzard days/cases
do i=1,ii
   do j=1,jj
      if(t2m(i,j,1,1)<-999)then
         cases_tot(i,j)=-999
         cases_annual(i,j)=-999
      else
         cases_tot(i,j)=0
         cases_annual(i,j)=0
         do y=1,yr
            cases(i,j,v)=0
         end do
      end if
   end do
end do
 
do i=1,ii
   do j=1,jj
      do yq=1,yr
         y=yq+yrbegin-1
         if(mod(y,4)==0.and.day_max==366)then
            if(v/=2100)then
               yrlength=366
            end if
         else
            yrlength=365
         end if
         if(day_max==360)then
            yrlength=day_max
         end if
         do day=1,yrlength
            if(t2m(i,j,day,yq)<lim_t2m.and.prec(i,j,day,yq)>lim_prec&
                 .and.gust(i,j,day,yq)>lim_gust)then
               cases(i,j,yq)=cases(i,j,yq)+1
            end if
         end do
      end do
   end do
end do
 
do i=1,ii
   do j=1,jj
      do y=1,yr
        cases_tot(i,j)=cases_tot(i,j)+cases(i,j,v)
        if(cases(i,j,v)>0)then
           cases_annual(i,j)=cases_annual(i,j)+1
        end if
     end do
  end do
end do
 
! Checking the results
do y=1,yr
write(6,*) v+yrbegin-1, cases(ipx,jpx,v)
end do
write(6,*) cases_tot(ipx,jpx),cases_annual(ipx,jpx)
 
! Saving the results
open(1,file='/RAIN/temporal/blizzards/&
bliz_EUR-44reg_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22E_v1_day_1971-2000.gr',&
     form='unformatted',status='unknown',access='direct',recl=4*ii*jj)
 
irec=1
do y=1,1
   write(1,rec=irec)((cases_tot(i,j),i=1,ii),j=1,jj)
   irec=irec+1
   write(1,rec=irec)((cases_annual(i,j),i=1,ii),j=1,jj)
   irec=irec+1
end do
close (1)
 
end program blizzards
