
global month_dict, totdays_dict
month_dict={'1':31,'2':28, '3':31, '4':30, '5':31, '6':30, '7':31, '8':31, '9':30, '10':31, '11':30, '12':31}
totdays_dict={'1':0,'2':31, '3':59, '4':90, '5':120, '6':151, '7':181, '8':212, '9':243, '10':273, '11':304, '12':334}

def fracitup(yr, mon, day, hr, mins):

	daystot = totdays_dict[str(int(mon))] + day + (hr + mins/60.) / 24.
	yrlen = 365
	if ((yr % 4 == 0) & (mon > 2)): daystot += 1
	if (yr % 4 == 0): yrlen = 366
	
	fracyr = daystot / yrlen
	return yr + fracyr
