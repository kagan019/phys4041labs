def mini(arr, lesser = lambda a,b : a < b):
    m = 0
    for i,x in enumerate(arr):
        if lesser(x,arr[m]):
            m = i
    return m
#a) 
def sort_asc(arr):
    newarr = []
    while len(arr) > 0:
        i = mini(arr)
        newarr += [arr[i]] 
        arr = arr[:i] + arr[i+1:]
    return newarr
#b)
def sort_asc_abs(arr):
    newarr = []
    while len(arr) > 0:
        i = mini(arr,lesser = lambda a, b : abs(a) < abs(b))
        newarr += [arr[i]] 
        arr = arr[:i] + arr[i+1:]
    return newarr

def main():
    testarr = [3,0,-34,543,-549,565,12,-3,-4938,-2,-968,533]
    print(str(sort_asc(testarr)))
    print(str(sort_asc_abs(testarr)))

if __name__ == "__main__":
    main()
