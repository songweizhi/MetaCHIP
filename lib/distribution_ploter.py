import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def plot_identity_list(identity_list, identity_cut_off, title, output_foler) :
    identity_list = sorted(identity_list)
    # get statistics
    match_number = len(identity_list)
    average_iden = float(np.average(identity_list))
    average_iden = float("{0:.2f}".format(average_iden))
    max_match = float(np.max(identity_list))
    min_match = float(np.min(identity_list))
    # get hist plot
    num_bins = 50
    plt.hist(identity_list,
             num_bins,
             alpha = 0.1,
             normed = 1,
             facecolor = 'blue')  # normed = 1 normalized to 1, that is probablity
    plt.title('Group: %s' % title)
    plt.xlabel('Identity')
    plt.ylabel('Probability')
    plt.subplots_adjust(left = 0.15)
    # get fit line
    density = gaussian_kde(identity_list)
    x_axis = np.linspace(min_match - 5, max_match + 5, 200)
    density.covariance_factor = lambda: 0.3
    density._compute_covariance()
    plt.plot(x_axis, density(x_axis))
    # add text
    x_min = plt.xlim()[0]  # get the x-axes minimum value
    x_max = plt.xlim()[1]  # get the x-axes maximum value
    y_min = plt.ylim()[0]  # get the y-axes minimum value
    y_max = plt.ylim()[1]  # get the y-axes maximum value
    # set text position
    text_x = x_min + (x_max - x_min)/5 * 3.8
    text_y_total = y_min + (y_max - y_min) / 5 * 4.4
    text_y_min = y_min + (y_max - y_min) / 5 * 4.1
    text_y_max = y_min + (y_max - y_min) / 5 * 3.8
    text_y_average = y_min + (y_max - y_min) / 5 * 3.5
    text_y_cutoff = y_min + (y_max - y_min) / 5 * 3.2
    # plot text
    plt.text(text_x, text_y_total, 'Total: %s' % match_number)
    plt.text(text_x, text_y_min, 'Min: %s' % min_match)
    plt.text(text_x, text_y_max, 'Max: %s' % max_match)
    plt.text(text_x, text_y_average, 'Mean: %s' % average_iden)
    plt.text(text_x, text_y_cutoff, 'Cutoff: %s' % identity_cut_off)
    if identity_cut_off != 'None':
        plt.annotate(' ',
                xy = (identity_cut_off, 0),
                xytext = (identity_cut_off, density(identity_cut_off)),
                arrowprops = dict(width = 0.5,
                                  headwidth = 0.5,
                                  facecolor = 'red',
                                  edgecolor = 'red',
                                  shrink = 0.02))
    else:
        pass
    # Get plot
    plt.savefig('%s/%s.png' % (output_foler, title), dpi = 300)
    plt.close()
