from solar_siblings_functions import *
import sys, getopt
from os.path import isfile

# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
# parse inputs from command line

process_bands = np.array([1])  # in range 1...4
read_ext = 0
process_obj_begin = 0
process_obj_end = -1
out_dir_suffix = ''

argv = sys.argv
if len(argv) > 1:
    # parse input options
    opts, args = getopt.getopt(argv[1:], '', ['bands=', 'ext=', 'obj_beg=', 'obj_end=', 'dir_suffix='])
    # set parameters, depending on user inputs
    print opts
    for o, a in opts:
        if o == '--bands':
            process_bands = np.array([np.int32(b) for b in a.split(',')])
            print 'Command line selected bands: ' + ','.join([str(pb) for pb in process_bands])
        if o == '--ext':
            read_ext = np.int8(a)
        if o == '--obj_beg':
            process_obj_begin = np.int64(a)
        if o == '--obj_end':
            process_obj_end = np.int64(a)
        if o == '--dir_suffix':
            out_dir_suffix = str(a)

d_wvl = 0.0
save_plots = True
output_differences = True  # so far available only for the first analysis step
min_wvl = min_wvl[process_bands-1]
max_wvl = max_wvl[process_bands-1]

GP_compute = True
save_gp_params = True
save_gp_params_read_append = True
n_threads = 15
w_m = 4
n_walkers = np.array([w_m*n_threads, w_m*n_threads, w_m*n_threads, w_m*n_threads])[process_bands-1]
n_steps = np.array([60, 60, 70, 80])[process_bands-1]

# evaluate spectrum
n_noise_samples = 750
noise_power = 0

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix_solar_ref = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'
# solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix, every_nth=8)

# data-table settings
data_date = '20180222'
galah_param_file = 'sobject_iraf_53_reduced_'+data_date+'.fits'
cannon_param_file = 'sobject_iraf_iDR2_180108_cannon.fits'

# select ok objects
print 'Reading and determining usable twilight flats for Solar statistics determination'
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist_newer.csv')
galah_param = Table.read(galah_data_input + galah_param_file)
cannon_param = Table.read(galah_data_input + cannon_param_file)
# join datasets and add some information to cannon parameters
cannon_param = join(cannon_param, galah_param['sobject_id','snr_c1_guess','snr_c2_guess','snr_c3_guess','snr_c4_guess'],
                    keys='sobject_id')

idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c2_guess'] > 200)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

# same for Cannon
idx_row_cannon = np.in1d(cannon_param['sobject_id'], galah_param[idx_rows]['sobject_id'])
idx_row_cannon = np.logical_and(idx_row_cannon, cannon_param['flag_cannon'] == 0)  # only unflagged flats for parameters

teff_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Teff_cannon'])
teff_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Teff_cannon'])
logg_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Logg_cannon'])
logg_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Logg_cannon'])
feh_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Feh_cannon'])
feh_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Feh_cannon'])
print 'Solar parameters - cannon:', teff_solar_c, '+/-', teff_solar_std_c, ',  ', logg_solar_c, '+/-', logg_solar_std_c, ',  ', feh_solar_c, '+/-', feh_solar_std_c

# manual parameter selection
idx_solar_like = (np.abs(cannon_param['Teff_cannon'] - teff_solar_c) <= 250) & \
                 (np.abs(cannon_param['Logg_cannon'] - logg_solar_c) <= 0.4) & \
                 (np.abs(cannon_param['Feh_cannon'] - feh_solar_c) <= 0.3)
# preform flag filtering if needed - later selection is currently implemented
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['flag_cannon'] >= 0)  # no flagging at this point
idx_solar_like = np.logical_and(idx_solar_like, np.bitwise_and(cannon_param['red_flag'], 64) == 0)  # only flats are taken out
# snr selection
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['snr_c2_guess'] >= 0)  # no snr limits
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['sobject_id'] > 140301000000000)  # leave out comissoning phase

n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

# -----------------------------------
# --------- Main program ------------
# -----------------------------------

solar_like_sobjects = cannon_param['sobject_id'][idx_solar_like]
sim_metrices = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'minkowski', 'wchebyshev', 'sqeuclidean', 'euclidean', 'chi2', 'EW', 'median_sep', 'sum', 'px_over', 'px_under']
sim_metrices_std = [m+'_std' for m in sim_metrices]
if GP_compute:
    sim_metrices_min = [m+'_min' for m in sim_metrices]
    sim_metrices_max = [m+'_max' for m in sim_metrices]
    sim_dtypes = ['float64' for i in range(4*len(sim_metrices))]
    sim_results = Table(names=np.hstack(('sobject_id', 'snr_spectrum', 'median_cont', 'median_cont_after', sim_metrices, sim_metrices_std, sim_metrices_min, sim_metrices_max)),
                        dtype=(np.hstack(('int64', 'float64', 'float64', 'float64', sim_dtypes))))
else:
    sim_dtypes = ['float64' for i in range(2*len(sim_metrices))]
    sim_results = Table(names=np.hstack(('sobject_id', 'snr_spectrum', 'median_cont', sim_metrices, sim_metrices_std)),
                        dtype=(np.hstack(('int64', 'float64', 'float64', sim_dtypes))))

bands_suffix = '_b'+''.join([str(b) for b in process_bands])
print bands_suffix
dir_suffix = '_p'+str(noise_power)+'_SNRsamples'+str(n_noise_samples)
dir_suffix += '_ext'+str(read_ext)
if OK_LINES_ONLY:
    dir_suffix += '_oklinesonly'
if not USE_SUBSAMPLE:
    dir_suffix += '_origsamp'
dir_suffix += out_dir_suffix
if GP_compute:
    move_to_dir(out_dir + 'Distances_Step2' + dir_suffix)
    gp_param_labels = ['amp_noise', 'rad_noise', 'amp_cont', 'rad_cont', 'cont_norm']
    txt_out_gp = 'GP_fit_res.txt'
    if isfile(txt_out_gp) and not save_gp_params_read_append:
        # do not rewrite original file
        txt = open(txt_out_gp, 'w')
        txt_out_header_str = 'sobject_id'
        for i_p_b in process_bands:
            txt_out_header_str += ',' + ','.join([p_l+'_b'+str(i_p_b) for p_l in gp_param_labels])
        txt.write(txt_out_header_str + '\n')
        txt.close()
else:
    move_to_dir(out_dir + 'Distances_Step1' + dir_suffix)

file_out_fits = 'solar_similarity'+bands_suffix+'.fits'
file_out_diff = 'solar_spectral_diff'+bands_suffix+'.csv'

# predetermined objects
solar_like_sobjects = [
140310003301028,140312003501087,140312003501360,140312004501092,140314004401272,140314004401277,140413003201328,140414004101334,140414004101367,140415002401342,140607001401278,140608001401244,140608002501303,140707001101347,140709003001184,140709003801149,140710000101230,140710000101252,140711001901181,140711001901389,140711002401078,140711002401243,140711002401316,140711002901069,140711002901153,140711002901268,140711003901161,140711003901328,140713001401273,140713001401308,140713004001112,140713004001188,140713004601016,140713004601346,140805002601115,140805003101138,140805003101303,140805003101338,140805003601343,140805004201051,140805004801129,140805004801137,140806000601014,140806000601206,140806002301027,140806002301392,140806002901279,140806003501105,140806004101159,140806004101186,140807000601003,140807004501336,140807005001072,140807005001277,140808002701286,140808002701299,140808002701338,140808002701394,140808003701104,140808003701261,140809002601201,140809003101032,140809003101322,140809003701037,140809003701146,140809003701186,140809004201211,140809004201296,140811003901103,140811003901267,140811004501035,140811004501343,140811005001245,140812003801322,140814003301151,140814003301261,141101001801006,141102002401182,141102002701072,141102002701267,141102002701342,141103003101379,141103003601164,141104004301324,141231004001331,150101004001347,150102003201119,150102003701187,150103002701017,150103003001048,150103003001128,150103004001229,150105002801352,150107004701146,150204002901346,150206004301295,150207002101314,150207003101166,150207004101059,150207004101353,150207005101103,150207005101176,150207005101334,150208003201286,150208003201339,150208004201342,150208005201146,150211004701179,150330002201348,150330002601242,150330002601306,150330002601320,150407001101120,150408004101169,150408004101278,150408004101283,150408005301076,150408005301093,150408005301116,150408005301184,150408005301259,150408005901078,150409002601317,150409003101129,150409003601089,150409004101294,150409005601390,150412002101135,150412002101208,150412002601212,150412002601280,150412005101226,150413003601248,150413003601344,150413003601355,150413005101096,150427002801068,150427004301259,150427004301273,150427004801068,150427004801275,150428000601336,150429001601268,150430002801394,150504003001047,150504003001111,150602002101246,150602002701341,150602003301072,150602003901271,150607005601279,150703003101181,150703003601249,150703004101375,150703005101366,150703005601061,150703005601062,150703005601082,150705001901052,150705005401225,150705005401272,150705005401363,150705005401364,150705005401377,150705006401129,150705006401314,150824002101315,150824002601134,150827005201080,150828002701298,150828003201089,150828003201374,150828003701040,150828004201034,150828004201212,150828004701096,150828005201075,150828005201385,150828005701069,150828005701331,150829002601283,150829003101225,150830002301065,150830002801072,150830004001040,150830004601175,150830005101021,150830005101144,150830005101337,150830005101391,150831003501011,150901000601082,151008003501085,151009001601351,151009001601363,151009003101282,151109002101254,151109003601038,151111002101059,151111002101116,151111002101236,151219001601078,151219002601083,151219002601248,151219003101397,151219003601105,151219003601245,151219003601298,151219003601330,151225002101224,151225002701158,151227005201142,151229004501012,151229004501170,151231002601074,151231002601098,160106002601047,160106003601037,160107001601085,160107002601011,160107002601186,160107003101157,160108002001207,160108003101163,160108003601152,160109003801085,160110002101345,160110003601039,160112001601056,160112002901151,160123002601062,160123002601079,160124003101351,160125001601112,160125002401126,160125003001258,160125004501038,160125004501256,160125004501389,160129004701062,160129004701355,160129005201070,160129005201118,160130003601121,160130004601023,160130005201082,160130005201356,160130005201357,160130006301234,160325002701048,160325002701090,160325002701275,160325003201071,160326001101111,160326001101395,160326002101077,160327003601167,160327003601210,160327003601277,160327003601386,160327004101056,160327004101248,160327004101343,160327004601025,160327004601337,160327006101355,160328000701292,160328003201282,160328004201051,160328004201352,160330002101179,160330002601095,160330002601230,160331002201390,160331002701393,160331004301386,160331004801091,160331004801301,160331004801396,160331005301074,160331005301201,160331005301347,160401002101259,160401003901022,160401003901032,160401003901215,160401004401051,160401004401097,160401004401121,160401004401123,160401004401168,160401005401170,160401005401203,160402004601106,160402005101147,160402005601084,160402005601358,160402005601383,160402006101178,160402006101388,160402006601248,160403003601392,160415001601014,160415002601198,160415002601342,160415003101244,160415003601021,160415004601223,160418003101331,160418004101386,160419002601052,160419002601161,160419003101371,160419003601012,160419003601116,160419003601127,160419003601201,160419004101106,160420003301342,160420004301031,160420004301188,160420004301224,160420004301292,160420005301119,160420006301085,160420006301234,160420006301373,160421002101006,160421002101145,160421003601190,160422002001351,160422002501096,160422002501382,160422003001267,160422003501066,160422003501107,160422004001270,160422004501162,160422004501312,160424002101194,160424002601134,160424003101290,160424004201264,160424004701369,160425001901071,160425001901089,160425001901273,160426004001348,160426004501264,160426004501385,160426004501395,160426005501042,160426006101034,160426006101240,160426006101290,160426006101338,160426006701393,160513001101030,160513001601022,160513002601086,160513003101265,160515002801245,160520002101297,160520002601048,160520002601397,160520004201020,160520004901083,160521004801353,160522002601049,160522003601193,160522004101119,160522004101313,160522006101314,160522006601193,160524002701116,160524002701328,160524004201306,160524004901098,160524005501254,160524005501270,160524006101090,160524006601258,160525002201033,160525002201197,160525002701071,160525003201154,160529001801017,160529003401069,160529005401085,160530002801021,160530003301365,160530003901026,160530003901054,160530005001359,160530005501077,160530005501111,160531001601046,160531001601307,160531004101086,160531004601362,160531005101256,160531005601362,160531006101153,160602001601307,160613001801082,160724003501032,160811002901141,160812002601392,160812003101017,160812003601094,160813002101147,160813002101295,160813002101314,160813003601074,160813003601151,160813003601290,160813004101092,160813004101136,160813004101239,160813004601107,160813005101077,160815002601022,160815003601239,160815004301182,160815004301313,160815004301316,160815004801066,160816002701321,160816003201280,160816003201323,160816003701030,160816003701288,160816004201351,160816004201365,160816004701036,160817002601023,160817002601053,160817002601285,160817003101217,160817003101276,160916001801076,160916001801111,160916001801263,160916003301120,160916003801373,160916004301036,160916004301152,160916004301242,160919001601330,160919003001290,160919003001321,160919003001365,160919003001367,160919004001328,160919005101057,160919005101083,160919005101297,160919005101351,160923005201179,161006003901387,161006004401018,161006004901009,161006004901099,161006004901385,161007002801316,161007002801397,161007003301394,161008002001043,161008002501018,161008003001018,161008003001092,161008003001222,161008003501226,161009002601018,161009002601361,161009003201386,161009003801064,161009003801163,161009004801144,161009004801344,161009005901171,161011003401069,161011003401169,161012002101063,161012002101185,161012002101227,161013002101025,161013002101244,161013002601207,161013004901384,161104002301129,161104002301257,161104002801361,161104002801394,161104003301091,161104003801323,161105003101019,161105003101345,161105004601015,161105004601374,161106003101031,161107001601102,161107001601132,161107001601133,161107004401071,161109003101231,161115002701019,161115002701213,161116001701213,161116001701314,161116002201127,161116002801375,161116003301345,161116003801308,161118002601069,161118002601176,161118002601376,161118003501203,161118004701023,161118004701221,161118004701364,161119002801292,161119003101279,161119003101362,161209001801389,161210004201073,161210004201315,161211003101098,161211003101387,161212002101397,161212002601031,161213003101374,161213003601035,161213004101187,161213004601295,161217002601138,161217003101387,161217004101075,161217004101234,161217004601271,161217006101046,161218002601077,161219003101122,161219004101199,161219004101351,161219004601026,161219005101228,161219005101286,170102001901072,170104002401006,170105003101089,170105003101138,170105004101367,170107004201309,170107004801107,170108003301041,170108004601041,170108004601140,170109001801030,170109001801166,170109002101081,170109002101106,170109002801065,170109003801047,170112002101294,170112002601348,170112003101027,170112003601284,170112003601298,170113001601045,170113001601220,170113001601289,170113002101159,170113002601194,170114001601154,170114002601280,170114004101216,170114004101329,170115001601273,170115003701382,170115004201141,170117003101044,170117003101267,170118002201216,170118003801289,170119002601393,170121002801292,170122002101017,170122002101113,170130002601157,170130003101184,170130003601019,170130004601325,170205003401264,170205003401299,170205003401309,170205003901001,170205003901143,170205004401235,170205004901019,170205004901178,170205004901192,170205004901388,170205005401120,170205005401129,170205005401143,170205005401275,170205005401391,170206003701346,170206004201399,170206004701139,170206005201202,170206005201287,170206005201289,170206005701018,170206005701110,170206005701345,170206005701369,170216003301325,170216003801347,170217001601332,170217002201016,170219001601323,170219001601353,170219002601038,170219002601324,170219003601351,170220002101322,170220004101045,170220004101215,170404001601387,170407002101137,170407002101264,170407002601374,170407003101281,170407003601111,170407003601187,170407003601304,170407004601072,170407004601205,170407005201023,170407005201138,170407005201381,170408004001026,170408004501048,170408004501192,170408004501227,170408004501368,170408005501159,170411002601201,170411003101025,170412002901210,170413002101206,170413002101361,170413002601080,170413002601218,170413003101029,170413003101069,170413004601076,170413005101082,170413005601032,170413005601051,170413005601334,170414002601089,170414003101035,170414003601014,170414003601069,170414004101039,170414005101057,170414005601224,170415002001080,170415002501325,170416003301103,170416003301117,170416003301342,170416003801070,170416003801092,170416004801107,170416004801120,170417003201330,170417005001332,170418001601351,170418002701078,170418002701261,170418003701097,170506002901026,170506002901071,170506002901139,170506003901066,170506004401274,170506006401014,170508001601184,170508002601001,170508004301053,170508004801312,170508005301075,170509002701171,170509002701386,170509004701096,170509005201063,170510003301024,170510005801062,170510006801361,170510007301226,170511000601263,170511001101015,170511001601237,170511004001234,170512000701245,170512001301074,170513004901118,170513004901180,170514002401099,170514003001180,170514003301001,170514003301011,170514003301312,170515002101186,170515003101035,170515003101036,170515003101127,170515003601030,170515006101166,170515006101289,170516000601373,170516001601139,170516002101273,170517001801023,170517002801389,170530001601155,170601003101179,170601003101218,170614002101078,170614002601314,170614003101390,170614004101220,170614004101305,170614004101381,170614004601030,170614004601045,170614004601055,170614004601061,170614005101328,170710002701354,170710003201082,170710003201145,170711002001243,170711002501066,170711002501169,170711003001134,170711004001131,170711004001148,170711004001339,170711004501033,170711005101185,170711005801034,170713004601335,170713004601346,170713005101388,170723003601118,170723005101385,170724003601147,170724004101366,170724004601055,170724005101385,170801002801345,170801004001010,170801004001139,170801004001215,170802002101151,170802003201111,170805003101077,170805003601318,170805003601381,170805005101295,170805005101392,170829001901039,170829003901119,170905001601158,170905002601072,170905003101147,170905003101295,170905003101357,170905003501094,170906002601322,170906003101307,170906003601056,170906003601147,170906003601159,170906003601210,170906003601295,170906003601357,170906004101067,170906004601391,170906005101282,170906005101299,170906005101318,170906005601108,170906005601141,170907002601114,170907002601378,170907003601280,170907004101051,170907004601162,170908001601149,170908002801274,170908002801276,170909001601114,170909002101136,170909002601291,170909002601340,170910001801313,170910002601041,170910003101074,170910003101081,170910003601009,170910003601321,170910004101188,170910004601141,170910005101082,170910005101268,170910005101316,170911002101145,170911002101388,170911002601049,170911003101334,170911003101383,170911004201047,170911004201366,170911004701065,170912002401132,170912002901303
]
print 'Number of pre-selected objects:', len(solar_like_sobjects)

# random subset objects from parameters selection
# n_rand = 25
# solar_like_sobjects = solar_like_sobjects[np.int64(np.random.rand(n_rand)*len(solar_like_sobjects))]

# first erase all results from previous processing runs
if output_differences:
    csv_diff = open(file_out_diff, 'w')
    csv_diff.close()

for s_obj in solar_like_sobjects[process_obj_begin:process_obj_end]:
    print 'Evaluating', s_obj
    galah_object = galah_param[galah_param['sobject_id'] == s_obj]
    # get spectra of all bands for observed objects
    # flux, wvl = get_spectra_dr52(str(s_obj), bands=[1, 2, 3, 4], root=dr52_dir, extension=read_ext, individual=False)
    flux, wvl, flux_std = get_spectra_dr52(str(s_obj), bands=process_bands, root=dr52_dir, extension=read_ext,
                                           individual=False, read_sigma=True)
    if len(flux) <= 0:
        continue
    if read_ext == 0:
        # normalize flux
        try:
            for i_c in range(len(process_bands)):
                # ------ NORM v1 - high order polynomial, many steps
                # flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=35, sigma_low=1.5, sigma_high=2.8, order=29,
                #                               n_min_perc=3.,  return_fit=False, func='poly')
                # ------ NORM v2 - the same as used in the process of reference Solar spectra construction
                norm_ok_mask = determine_norm_mask(wvl[i_c], norm_bad_ranges)
                flux[i_c] = spectra_normalize(wvl[i_c] - np.mean(wvl[i_c]), flux[i_c], fit_mask=norm_ok_mask,
                                              steps=15, sigma_low=2., sigma_high=3., order=5, n_min_perc=5.,
                                              return_fit=False, func='cheb')
                # # additional normalization step with symmetric sigma rejection intervals to cancel out noise
                # flux[i_c] = spectra_normalize(wvl[i_c] - np.mean(wvl[i_c]), flux[i_c], fit_mask=norm_ok_mask,
                #                               steps=15, sigma_low=2.5, sigma_high=2.5, order=1, n_min_perc=5.,
                #                               return_fit=False, func='poly')
            # apply computed rv shift to the spectrum
            rv_shift = galah_object['rv_guess_shift']
            wvl *= (1 - rv_shift / 299792.458)
        except:
            print ' -> Something wrong with spectra or reading'
            continue

    # compute guess like snr for particular spectrum and observed region
    # get absorption features indices
    idx_lines_mask = get_linelist_mask(np.hstack(wvl))
    wvl_all_abs = np.hstack(wvl)#[idx_lines_mask]
    flx_all_abs = np.hstack(flux)#[idx_lines_mask]
    # median signal at selected abundance wavelength pixels
    snr_signal = np.nanmedian(flx_all_abs)
    # determine actual snr of generated noise at selected pixels - guess like
    snr_noise = 1.4826 / np.sqrt(2) * np.nanmedian(np.abs(flx_all_abs[1:] - flx_all_abs[:-1]))
    snr_guesslike = snr_signal / snr_noise
    # print 'SNRs:', galah_object['snr_c' + str(i_c + 1) + '_guess'].data[0], snr_guesslike

    # determine continuum-like pixels
    min_cont_level = 0.98
    solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=4)
    idx_cont_px = solar_flx > min_cont_level
    flx_all_abs_res = spectra_resample(flx_all_abs, wvl_all_abs, solar_wvl, k=1)
    idx_cont_px = np.logical_and(idx_cont_px, np.isfinite(flx_all_abs_res))
    cont_median_solar = np.nanmedian(solar_flx[idx_cont_px])
    cont_median_flx = np.nanmedian(flx_all_abs_res[idx_cont_px])
    cont_median_dif = cont_median_solar - cont_median_flx

    pix_ref = list([])
    pix_ref_noise = list([])
    pix_ref_cont = list([])  # continuum pixels for continuum offset determination
    pix_spec = list([])
    pix_std = list([])
    if GP_compute:
        gp_final_res = list([])
        # Start GP process for every band in spectrum independently

        for i_c in range(len(process_bands)):
            evaluate_band = process_bands[i_c]
            # first prepare reference data
            solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref,
                                                  every_nth=every_nth_solar_pixel[evaluate_band - 1])
            # band wvl mask
            idx_ref = get_band_mask(solar_wvl, evaluate_band)
            # generate mask of pixels used in comparison
            idx_lines_mask = get_linelist_mask(solar_wvl)

            # define subset of spectra to be compared to reference solar spectrum
            abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]

            flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
            flux_std_b_res = spectra_resample(flux_std[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)

            # correct flux values if needed
            flux_b_res[flux_b_res > 1.2] = 1.2
            flux_b_res[flux_b_res < 0] = 0.

            # determine spectrum difference and its variance
            diff = (solar_flx[idx_ref] - flux_b_res)
            diff_var = np.nanvar(diff)

            skip_gp_computation = False
            # check if GP results are already available:
            if isfile(txt_out_gp) and save_gp_params_read_append:
                try:
                    gp_precom_fit = Table.read(txt_out_gp, format='ascii.csv')
                    idx_line = np.where(gp_precom_fit['sobject_id'] == s_obj)[0]
                    if len(idx_line) == 1:
                        # we found a match, read it for current band
                        skip_gp_computation = True
                        gp_res_read_cols = [gp_p_l + '_b' + str(evaluate_band) for gp_p_l in gp_param_labels]
                        kernel_fit = gp_precom_fit[idx_line][gp_res_read_cols].to_pandas().values[0]
                        print ' GP emcee parameters restored'
                except:
                    print ' Problem restoring GP emcee parameters'

            if not skip_gp_computation:
                # determine kernel parameters trough emcee fit
                print ' Running emcee'
                # emcee_fit_px = 100
                rad_noise_init = [0.0018, 0.0030, 0.0040, 0.0055][evaluate_band - 1]  # start process with different initial values for every band
                sampler, fit_res, fit_prob = fit_gp_kernel([diff_var/2., rad_noise_init, 1e-5, 16, 1.],
                                                           solar_flx[idx_ref], flux_b_res, solar_wvl[idx_ref],
                                                           data_std=None,
                                                           # data_std=flux_std_b_res,  # takes even longer to compute GP
                                                           nwalkers=n_walkers[i_c], n_threds=n_threads, n_burn=n_steps[i_c],
                                                           exit_lnp=10, normal_dist_guess=False)
                # walker prob plot
                if save_plots:
                    print(" Plotting walker probabilities")
                    walkers_prob = sampler.lnprobability/len(flux_b_res)
                    for i_w in range(walkers_prob.shape[0]):
                        plt.plot(walkers_prob[i_w, :], lw=0.3)
                    walkers_prob = walkers_prob.flatten()                   # without this correction numpy
                    walkers_prob = walkers_prob[np.isfinite(walkers_prob)]  # percentile may return incorrect -inf value
                    plt.ylim((np.percentile(walkers_prob, 1), np.percentile(walkers_prob, 99.9)))
                    plt.savefig(str(s_obj) + '_gp-lnprob_b' + str(i_c + 1) + '.png', dpi=250)
                    # plt.show()
                    plt.close()

                sampler_chain_vals = sampler.flatchain
                kernel_fit = np.median(sampler_chain_vals, axis=0)  # flatchain holds parameters of all emcee steps

                # corner plot of parameters
                if save_plots:
                    c_fig = corner.corner(sampler.flatchain, truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
                                          labels=gp_param_labels, bins=30)
                    c_fig.savefig(str(s_obj) + '_corner_b' + str(i_c + 1) + '.png', dpi=200)
                    plt.close(c_fig)

            # add fitted values to the resulting table
            gp_final_res.append(kernel_fit)

            # create a gaussian process that will be used for the whole spectra
            gp = george.GP(get_kernel(kernel_fit[:-1]))
            gp.compute(solar_wvl[idx_ref])
            gp_noise_pred = gp.sample(size=n_noise_samples)
            gp_noise_pred += spectrum_offset_norm(kernel_fit[-1:], solar_flx[idx_ref]) - solar_flx[idx_ref]

            if save_plots:
                plt.plot(solar_flx[idx_ref], c='red', lw=0.5)
                for i_pred in range(20):
                    plt.plot(solar_flx[idx_ref] + gp_noise_pred[i_pred, :], c='black', alpha=0.15, lw=0.3)
                plt.plot(flux_b_res, c='blue', lw=0.5)
                plt.ylim((0.4, 1.1))
                # plt.show()
                plt.savefig(str(s_obj)+'_gp_b'+str(i_c+1)+'.png', dpi=400)
                plt.close()

            # determine continuum-like pixels after GP fitting was performed
            pix_ref_cont.append((solar_flx[idx_ref] + gp_noise_pred)[:, np.where(solar_flx[idx_ref] > min_cont_level)[0]])

            # store results for current band
            pix_ref.append(solar_flx[idx_ref][abs_lines_cols])
            pix_ref_noise.append(gp_noise_pred[:, abs_lines_cols])
            pix_spec.append(flux_b_res[abs_lines_cols])
            pix_std.append(flux_std_b_res[abs_lines_cols])

        # determine continuum median difference after GP fitting was performed
        cont_median_dif_after = np.median(np.hstack(pix_ref_cont)) - cont_median_flx

        # save fit res if they are not already in it
        if save_gp_params and not skip_gp_computation:
            txt = open(txt_out_gp, 'a')
            gp_res_string = str(s_obj) + ',' + ','.join([str(v) for v in np.array(gp_final_res).flatten()])
            txt.write(gp_res_string + '\n')
            txt.close()

    else:
        for i_c in range(len(process_bands)):
            evaluate_band = process_bands[i_c]
            # first prepare reference data
            solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref,
                                                  every_nth=every_nth_solar_pixel[evaluate_band - 1])
            # band wvl mask
            idx_ref = get_band_mask(solar_wvl, evaluate_band)
            # generate mask of pixels used in comparison
            idx_lines_mask = get_linelist_mask(solar_wvl)

            # define subset of spectra to be compared to reference solar spectrum
            abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]
            # print flux[i_c], wvl[i_c], solar_wvl[idx_ref]
            flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
            flux_std_b_res = spectra_resample(flux_std[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)

            if n_noise_samples <= 0:
                # do not add any noise to used reference spectra
                snr_noise_pred = np.zeros((1, len(flux_b_res)))
            else:
                # generate poissonian noise to make a spectrum with snr into a spectrum with target snr
                snr_ref = np.inf
                snr_spectrum = galah_object['snr_c' + str(i_c + 1) + '_guess'].data
                snr_sigma = np.sqrt((1.0 / snr_spectrum) ** 2)  # - (1.0 / snr_ref) ** 2)
                snr_noise_pred = np.random.poisson((1.0 / snr_sigma)**2, size=(n_noise_samples, len(flux_b_res)))
                snr_noise_pred = snr_noise_pred / ((1.0 / snr_sigma)**2) - 1.

            # store results for current band
            pix_ref.append(solar_flx[idx_ref][abs_lines_cols])
            pix_spec.append(flux_b_res[abs_lines_cols])
            pix_std.append(flux_std_b_res[abs_lines_cols])
            pix_ref_noise.append(snr_noise_pred[:, abs_lines_cols])

    # compute different distance measurements
    pix_ref = np.hstack(pix_ref)
    pix_ref_noise = np.hstack(pix_ref_noise)
    pix_spec = np.hstack(pix_spec)
    pix_std = np.hstack(pix_std)

    if not evaluate_spectrum(pix_spec, flux_std):
        continue

    # iterate and add noise to observed spectrum
    n_distances_compute = np.max([n_noise_samples, 1])
    spectrum_distances = np.zeros((n_distances_compute, len(sim_metrices)))
    if save_plots:
        plt.figure(1, figsize=(12, 7))
        axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
        axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])
    for i_snr in range(n_distances_compute):
        # determine weights for the distance computation (different approaches)
        spectrum_distances[i_snr, :] = compute_distances(pix_spec, pix_std, pix_ref_noise[i_snr, :] + pix_ref, d=noise_power)
        if save_plots:
            axSpectra.plot(pix_ref_noise[i_snr, :] + pix_ref, lw=0.2, alpha=0.01, c='blue')
        if output_differences and not GP_compute:
            csv_diff = open(file_out_diff, 'a')
            if os.path.getsize(file_out_diff) == 0:  # size of file is zero -> add wavelength header info
                csv_diff.write('0,'+','.join([str(sw) for sw in solar_wvl[idx_ref][abs_lines_cols]])+'\n')
            diff_csv_string = ','.join([str(pf) for pf in (pix_spec - (pix_ref_noise[i_snr, :] + pix_ref))])
            csv_diff.write(str(s_obj)+','+diff_csv_string+'\n')
            csv_diff.close()

    # add agregated results to final table
    if GP_compute:
        sim_results.add_row(np.hstack([s_obj, snr_guesslike, cont_median_dif, cont_median_dif_after, np.nanmean(spectrum_distances, axis=0),
                            np.nanstd(spectrum_distances, axis=0), np.nanmin(spectrum_distances, axis=0), np.nanmax(spectrum_distances, axis=0)]))
    else:
        sim_results.add_row(np.hstack([s_obj, snr_guesslike, cont_median_dif, np.nanmean(spectrum_distances, axis=0), 
                            np.nanstd(spectrum_distances, axis=0)]))

    if save_plots:
        axSpectra.plot(pix_ref, c='black', lw=0.5)
        axSpectra.plot(pix_spec, c='blue', lw=0.5)
        axDiff.axhline(y=0, c='black', lw=0.5)
        axDiff.plot(pix_ref-pix_spec, c='blue', lw=0.5)
        axSpectra.set(ylim=(0.3, 1.15))
        axDiff.set(ylim=(-0.05, 0.05))
        axSpectra.set_title('SNR from abs lines: {:.2f}'.format(snr_guesslike))
        plt.savefig(str(s_obj) + '_' + str(galah_object['snr_c2_guess'].data[0])+bands_suffix+'.png', dpi=300)
        plt.close()

# check output file with results
if os.path.isfile(file_out_fits):
    os.remove(file_out_fits)
sim_results.write(file_out_fits)

'''
print sim_results5
print ''
sobj_id_like = sim_results[np.argsort(sim_results['chi2'])[:75]]['sobject_id']
print ','.join([str(s) for s in sobj_id_like])

print ''
sobj_id_dislike = sim_results[np.argsort(sim_results['chi2'])[-75:]]['sobject_id']
print ','.join([str(s) for s in sobj_id_dislike])

# output a plot of the most solar like spectra
for i_b in range(1, 5):
    for s_obj in sobj_id_like:
        flux, wvl = get_spectra_dr52(str(s_obj), bands=[i_b], root=dr52_dir, extension=read_ext)
        if read_ext == 0:
            # apply the same normalization as in the process of creation of master solar spectrum
            flux[0] = spectra_normalize(wvl[0], flux[0], steps=35, sigma_low=1.5, sigma_high=2.8, order=29,
                                          n_min_perc=3., return_fit=False, func='poly')
            # apply computed rv shift to the spectrum
            rv_shift = galah_param[galah_param['sobject_id'] == s_obj]['rv_guess_shift']
            wvl *= (1 - rv_shift / 299792.458)
        plt.plot(wvl[0], flux[0], lw=0.2, c='blue', alpha=0.02)
    plt.plot(solar_wvl, solar_flx, lw=0.2, c='black')
    plt.xlim((min_wvl[i_b-1], max_wvl[i_b-1]))
    plt.ylim((0.4, 1.1))
    plt.savefig('similar_spectra_b'+str(i_b)+'.png', dpi=1000)
    plt.close()
'''
