import numpy as np

def fun_practice(arr_of_reals):
    return np.sign(sum(arr_of_reals))

def main():
    s = fun_practice(np.array([-1,6,-34,-60,32]))
    print(s)

if __name__ == "__main__":
    main()
