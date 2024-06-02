import seaborn as sns
import matplotlib.pyplot as plt


def plot_heatmap(mydf, save="no", myoutput=None, vmin = 3, vmax = 20, cmap = None, title =None):
     c = (save.lower() == "no") or (save.lower() == "yes")
     # print(c)
     assert c, "ERROR: save must be one of yes or no"

     fig, ax = plt.subplots(figsize=(24,18))
     sns.heatmap(mydf, vmin = vmin, vmax = vmax, cmap = cmap)

     y1_1 = list(mydf.index.values).index("26") + 0.5
     y1_2 = list(mydf.index.values).index("39")+ 0.5
     x1_1 = list(mydf.columns.values).index("26") + 0.5
     x1_2 = list(mydf.columns.values).index("39") + 0.5

     y2_1 = list(mydf.index).index("55") + 0.5
     y2_2 = list(mydf.index).index("66") + 0.5
     x2_1 = list(mydf.columns).index("55") + 0.5
     x2_2 = list(mydf.columns).index("66") + 0.5

     y3_1 = list(mydf.index).index("104") + 0.5
     y3_2 = list(mydf.index).index("118") + 0.5
     x3_1 = list(mydf.columns).index("104") + 0.5
     x3_2 = list(mydf.columns).index("118") + 0.5

     y44_1 = list(mydf.index).index("43") + 0.5
     y44_2 = list(mydf.index).index("45") + 0.5
     x44_1 = list(mydf.columns).index("43") + 0.5
     x44_2 = list(mydf.columns).index("45") + 0.5

     poly_coords_cdr11 = [(x1_1, y1_1), (x1_1, y1_2),(x1_2, y1_2), (x1_2, y1_1)]
     poly_coords_cdr12 = [(x1_1, y2_1), (x1_1, y2_2),(x1_2, y2_2), (x1_2, y2_1)]
     poly_coords_cdr13 = [(x1_1, y3_1), (x1_1, y3_2),(x1_2, y3_2), (x1_2, y3_1)]
     poly_coords_cdr21 = [(x2_1, y1_1), (x2_1, y1_2),(x2_2, y1_2), (x2_2, y1_1)]
     poly_coords_cdr22 = [(x2_1, y2_1), (x2_1, y2_2),(x2_2, y2_2), (x2_2, y2_1)]
     poly_coords_cdr23 = [(x2_1, y3_1), (x2_1, y3_2),(x2_2, y3_2), (x2_2, y3_1)]
     poly_coords_cdr31 = [(x3_1, y1_1), (x3_1, y1_2),(x3_2, y1_2), (x3_2, y1_1)]
     poly_coords_cdr32 = [(x3_1, y2_1), (x3_1, y2_2),(x3_2, y2_2), (x3_2, y2_1)]
     poly_coords_cdr33 = [(x3_1, y3_1), (x3_1, y3_2),(x3_2, y3_2), (x3_2, y3_1)]

     poly_coords_Q44 = [(x44_1, y44_1), (x44_1, y44_2),(x44_2, y44_2), (x44_2, y44_1)]

     ax.add_patch(
          plt.Polygon(poly_coords_cdr11, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr12, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr13, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr21, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr22, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr23, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr31, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr32, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )
     ax.add_patch(
          plt.Polygon(poly_coords_cdr33, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )

     ax.add_patch(
          plt.Polygon(poly_coords_Q44, color ='black',
                    fill = None, lw = 2,
                    alpha = 1) )

     plt.xticks(rotation=90)
     if title != None:
        plt.title(title)
     if save == "yes":
          assert myoutput, "ERROR: please assign output folder and name"
          plt.savefig(myoutput, format = "svg")
     plt.show()