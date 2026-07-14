from pwr_f2_test import pwr_f2_test
import matplotlib.pyplot as plt
import numpy as np

outdir = '../figures/'

R2s = np.array([.1,.2,.3,.4,.5,.6,.7,.8,.9])
print(R2s)

f2s = R2s / (1-R2s)
f2s = [np.round(f2,3) for f2 in f2s]
print(f2s)

num_obs = [pwr_f2_test(u = 2, f2 = f2, sig_level = 0.05, power = 0.8)[-1] for f2 in f2s]
N = 13

print(f2s)
print(num_obs)

ax = plt.subplot()
ax.plot(num_obs, R2s, c = 'tab:blue', label = r'$R^2$')
ax.set_ylabel(r'Detectable $R^2$', c = 'tab:blue')
ax2 = ax.twinx()
ax2.plot(num_obs, f2s, c = 'tab:orange', label = r'$f^2$')
ax2.set_ylabel(r'Detectable $f^2$', c = 'tab:orange', rotation = 270)
ax2.axhline(.02, c = 'grey', ls = ':')
ax2.axhline(.15, c = 'grey', ls = ':')
ax2.axhline(.35, c = 'grey', ls = ':')
ax2.text(93, .02, 'Small effect size', c='grey', fontsize=10, ha = 'right')
ax2.text(93, .15, 'Medium effect size', c='grey', fontsize=10, ha = 'right')
ax2.text(93, .35, 'Large effect size', c='grey', fontsize=10, ha = 'right')
ax2.set_yscale('log')
ax.set_xlabel('Number of observation')
plt.axvline(N, c = 'r', ls = '--')
ax.text(N+1,.9, 'N = ' + str(N), c = 'r')
plt.title(r'80% power, $\alpha=0.05$')
l1,h1 = ax.get_legend_handles_labels()
l2,h2 = ax2.get_legend_handles_labels()
plt.legend(l1+l2, h1+h2)
plt.savefig(outdir + 'Power_quadratic_13epitopes.png')
plt.show()

num_obs = [pwr_f2_test(u = 1, f2 = f2, sig_level = 0.05, power = 0.8)[-1] for f2 in f2s]
N = 22

ax = plt.subplot()
ax.plot(num_obs, R2s, c = 'tab:blue', label = r'$R^2$')
ax.set_ylabel(r'Detectable $R^2$', c = 'tab:blue')
ax2 = ax.twinx()
ax2.plot(num_obs, f2s, c = 'tab:orange', label = r'$f^2$')
ax2.set_ylabel(r'Detectable $f^2$', c = 'tab:orange', rotation = 270)
ax2.axhline(.02, c = 'grey', ls = ':')
ax2.axhline(.15, c = 'grey', ls = ':')
ax2.axhline(.35, c = 'grey', ls = ':')
ax2.text(75, .02, 'Small effect size', c='grey', fontsize=10, ha = 'right')
ax2.text(75, .15, 'Medium effect size', c='grey', fontsize=10, ha = 'right')
ax2.text(75, .35, 'Large effect size', c='grey', fontsize=10, ha = 'right')
ax2.set_yscale('log')
ax.set_xlabel('Number of observation')
plt.axvline(N, c = 'r', ls = '--')
ax.text(N+1,.9, 'N = ' + str(N), c = 'r')
plt.title(r'80% power, $\alpha=0.05$')
l1,h1 = ax.get_legend_handles_labels()
l2,h2 = ax2.get_legend_handles_labels()
plt.legend(l1+l2, h1+h2)
plt.savefig(outdir + 'Power_linear_22epitopes.png')
plt.show()