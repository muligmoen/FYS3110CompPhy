import matplotlib.pyplot as plt

exact = 0.192765711

N = [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000]

MCC = [0.000000008, 0.001133467, 0.024285372, 0.014531879,
       0.194244256, 0.117016553, 0.206257811, 0.185448419]
MCCsigma = [0.000000007, 0.001125540, 0.023123862, 0.005365095,
            0.086522718, 0.027317216, 0.040739848, 0.012420672]
  
MCP = [0.000376852, 0.231511633, 0.157525515, 0.212161377,
       0.197445741, 0.206297434, 0.198070682, 0.200285570]
MCPsigma = [0.000145215, 0.126265106, 0.027329668, 0.016812296,
            0.004107642, 0.007374732, 0.000595095, 0.000732391]
  
MCI = [0.240731444, 0.412266740, 0.169150179, 0.190142914, 
       0.193869185, 0.193403734, 0.195505813, 0.195270296]
MCIsigma = [0.133228794, 0.177025271, 0.025264961, 0.006359748,
            0.002517490, 0.000843820, 0.000325820, 0.000270642]
  
  


ax = plt.subplot(111)

ax.set_xscale("log")

plt.plot([N[0], N[-1]], [exact, exact], linewidth=2, label='Exact')
plt.errorbar(N, MCC, yerr=MCCsigma, linewidth=2, label='MC cartesian')
plt.errorbar(N, MCP, yerr=MCPsigma, linewidth=2, label='MC polar')
plt.errorbar(N, MCI, yerr=MCCsigma, linewidth=2, label='MC importance')
#plt.plot(N, rel_err[1], linewidth=2, label='MC polar')
#plt.plot(N, rel_err[2], linewidth=2, label='MC importance')


plt.xlabel('N')
plt.ylabel('Integral estimate')
plt.xlim((N[0], N[-1]))
plt.ylim((0, 0.5))
plt.legend()

plt.savefig('MC_error.png',dpi=400, bbox_inches='tight')

plt.show()

