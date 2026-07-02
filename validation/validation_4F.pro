;;@package phosim
;; @file validation_temp.pro
;; @brief validation 4F
;;
;; @brief Created by:
;; @author John R. Peterson (Purdue)
;;
;; @brief Modified by:
;; @author Caleb Remocaldo (Purdue)
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.

;; Should look like a parabola roughly
;; Photodiode, canon, yoachim, all junk calibrations

pro validation_4f, nnn, vers, value, tolerance_low, tolerance_high, task, name, unit, comparison
print, 'Task 4F'
  
!p.multi = [0, 2, 2]
aa = findgen(33)/32.*2*!Pi
usersym, cos(aa), sin(aa), /fill

spawn, '/bin/ls wiyn_odi_e_* > wiynFiles.txt'
readcol, 'realData.txt', mjd, intensity, filter, format = '(D), (D), (A)', delimiter=' ', /silent
readcol, 'wiynFiles.txt', wiynFiles, format = '(A)', /silent

;'w' at the beginning of a variable signifies data from the simulated wiyn images
wmjd = dblarr(N_elements(wiynFiles))
wintensity = dblarr(N_elements(wiynFiles))
wfilter = strarr(N_elements(wiynFiles))
wfilt = strarr(N_elements(wiynFiles))

for i = 0L, N_elements(wiynFiles)-1 do begin
  print, "Reading in data from " + wiynFiles[i]
  
  data = mrdfits(wiynFiles[i], /silent)
  fitsHeader = headfits(wiynFiles[i], /silent)

  wfilter(i) = sxpar(fitsHeader, 'FILTER')
  wmjd[i] = sxpar(fitsHeader, 'MJD-OBS')
  wexptime = sxpar(fitsHeader, 'EXPTIME')

  ;should be counts per second per pixel
  wintensity[i] = ((mean(data(70:3950, 30:3610))*100.0/wexptime)>1.0)*wexptime
endfor

minMjd = floor(min(mjd))

print, "Determining nights"
;Determine nights
sortedMjd = sort(mjd)
pas = mjd[sortedMjd]
nightArr = pas[0]
night = 1L
for i = 0L, N_elements(pas)-2 do begin
   if pas(i+1)-pas(i) gt 0.25 then begin
      nightArr = [nightArr, pas(i+1)]
      if night eq 1 then dayArr = pas[i] else dayArr = [dayArr, pas[i]]
      night = night+1L
   endif
endfor

;If there is only one day, use the max MJD value
if N_elements(dayArr) eq 0 then begin
  dayArr = [max(mjd)]
endif else begin
  dayArr = [dayArr, max(pas)]
endelse

!p.multi = [0, 1, 2]

print, "Plotting data..."
filt = ['u', 'g', 'i', 'r', 'z']
colour = [50, 150, 200, 230, 0] ; Blue, green, orange, red, black

totalPoints = 0
wtotalPoints = 0
for ii = 0, N_elements(dayArr)-1 do begin
    ;Adjust x-values to be delta MJD, not just MJD
    shiftedNight = nightArr[ii]-minMjd
    shiftedDay = dayArr[ii]-minMjd

    ;Find all places where the image's shifted MJD value is within the current night cycle
    g = where(mjd-minMjd ge shiftedNight and mjd-minMjd le shiftedDay)
    wg = where(wmjd-minMjd ge shiftedNight and wmjd-minMjd le shiftedDay)
    
    realminy = floor(min([intensity(g)])/2.5)
    realmaxy = ceil(max([intensity(g)])*2.5)
    wminy = floor(min([wintensity(wg)])/2.5)
    wmaxy = ceil(max([wintensity(wg)])*2.5)

    miny = wminy<realminy
    maxy = wmaxy>realmaxy

  plot, [shiftedNight-.025, shiftedDay+.025], [miny, maxy], /nodata, /ylog, xtitle = '!4D!3MJD' + ' from ' + string(minMjd), ytitle = 'Counts/exp/pixel', /xstyle, /ystyle ;may need /ylog
  
  oplot, mjd[g]-minMJD, intensity[g], color = 50, psym = 5
  ;oplot, wmjd[wg]-minMJD, wintensity[wg], color = 200, psym = 4

  ;Plot real points in order of filter, checking for empty arrays
  for j = 0L, N_elements(filt)-1 do begin
     temp = where(filter[g] eq filt[j])
     if temp[0] ne -1 then begin
      oplot, mjd[g[temp]]-minMjd, intensity[g[temp]], color=colour[j], psym = 5 ;Triangles
      totalPoints += N_elements(temp)
     endif else begin
      print, "Plot " + string(ii+1) + ": no real data for filter " + string(filt[j])
     endelse
  endfor

  ;Plot simulated points in order of filter, checking for empty arrays
  for x = 0L, N_elements(filt)-1 do begin
     wtemp = where(wfilter[wg].contains(filt[x]))
     if wtemp[0] ne -1 then begin
      oplot, wmjd[wg[wtemp]]-minMjd, wintensity[wg[wtemp]], color=colour[x], psym = 4 ;Diamonds
      wtotalPoints += N_elements(wtemp)
     endif else begin
      print, "Plot " + string(ii+1) + ": no phosim data for filter " + string(filt[x])
     endelse
  endfor
  
  ;Only make legend on first plot
  if ii eq 0 then begin
    items = ["Observation data", "PhoSim data"]
    sym = [5, 4]
    ;For some reason, IDL does not automatically compile astron's 'legend' procedure.  To fix this, just compile it manually and the problem should go away.
    legend, items, psym=sym, /right
  endif
endfor
print, "Plotted " + string(totalPoints) + " real data points and " + string(wtotalPoints) + " simulated points."
END
