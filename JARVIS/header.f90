program CMOL
implicit double precision (a-h,o-z)
character(10) date,ztime
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'   OOOO    OO    OOOOOO  OOO OOO  OOOOO   OOOOO'
write(761,*)'     O      O     O    O  O   O     O    O     O'
write(761,*)'     O      O     O    O  O   O     O    O'
write(761,*)'     O     O O    O    O  O   O     O    O'
write(761,*)'     O     O O    OOOOO    O O      O     OOOOO'
write(761,*)'     O    O   O   O  O     O O      O          O'
write(761,*)' O   O    OOOOO   O  O     O O      O          O'
write(761,*)' O   O    O   O   O   O     O       O    O     O'
write(761,*)'  OOO    OOO OOO OOO  OO    O     OOOOO   OOOOO'
write(761,*)''
write(761,*)'          Product of Stark Industries          '
write(761,*)''
write(761,*)'    -----------------------------------------------'
write(761,*)'    |      WRITTEN BY DeCARLOS TAYLOR (March2009) |'
write(761,*)'    |            Revised April 12 (2009)          |'
write(761,*)'    |            Revised February (2012)          |'
write(761,*)'     -----------------------------------------------'
write(761,*)''
call date_and_time(date,ztime)
write(761,*)'Calculation started at: ',ztime(1:2),':',ztime(3:4),':',ztime(5:6) &
,' on ',date(5:6) ,'/', date(7:8),'/',date(1:4)
call cpusec(time1)
close(761)
call arl_input
call cpusec(time2)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)''
write(761,*)''

write(761,200)time2-time1

200 format(' Total job time (sec) =  ',F10.3)
write(761,*)''
call date_and_time(date,ztime)
write(761,*)'Calculation completed at: ',ztime(1:2),':',ztime(3:4),':',ztime(5:6) &
,' on ',date(5:6) ,'/', date(7:8),'/',date(1:4)
call space(1)
write(761,*)' Thanks for running DT-POLY :-)'
close(761)


end program cmol
