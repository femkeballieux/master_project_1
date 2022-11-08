All the code here is kind of a mess, so here is a short description of what everything does:

Code I inherited from Martje:
1. Martje_data_selec_isolated: takes in the LoTSS survey and checks whether the sources are isolated in here
2. Martje_data_selec_resolved_sourcetype: Takes in the isolated LoTSS sources and removes all those that have a weird sourcetype and are resolved
3. Martje_crossmatch_radio_surveys: crossmatch the resolved isolated LoTSS with the other radio surveys
4. Martje_make_tab_spectrum: cleans up the table from the crossmatching such that for all the surveys we are only left with flux, RA, Dec, fluxerror. 
            Also makes SEDs if you want it to
5. seds_plot_func_extreme: a module that contains functions for above
6. Martje_plot_colour_colour: makes a colour-colour plot from the table as given by nr. 4

Code I wrote myself:
7. change_columnnames: changes the columnnames because I cant type
8. crossmatch_master_inband_LoLSS: crossmatches the master sample we had above with also the LoLSS inband spectra
9. make_tab_master_inband_LoLSS: cleans up the table from above to get the right column names and units
10. make_plot_singular: is able to make a plot for a single, or for all PS sources from the table above

Because LoLSS is giving me a headache
11. crossmatch_old_new_master
12. compare_old_new_master 
13. crossmatch_old_new_LoLSS
14. compare_old_new_LoLSS

Now that we are done with LoLLS
15. make_tab_PS: so we can work with only the PS sources for example for VLASS things

others:
gpscssmodels: From Joe, be careful because there is a small mistake in there
