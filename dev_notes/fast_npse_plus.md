## Make NPSE+ fast again

Single threaded, Newton-Rhapson with matrices.
(5x5) pavement
``` 
[ Info: _________________________________________________
[ Info: Pavement time    = 1829.976
[ Info: % time in solver = 1664.397, 91 % of pavement time
[ Info: Single tile time = 66.576
[ Info: _________________________________________________
```

Not using collapse shortcut
(10x10)
```
[ Info: _________________________________________________
[ Info: Pavement time    = 3093.837
[ Info: % time in solver = 2853.452, 92 % of pavement time
[ Info: Single tile time = 28.535
[ Info: _________________________________________________
```


(5x5) multithreaded
```
[ Info: _________________________________________________
[ Info: Pavement time    = 1154.177
[ Info: % time in solver = 1111.748, 96 % of pavement time
[ Info: Single tile time = 44.470
[ Info: _________________________________________________
```

(20x20) multi
```
[ Info: _________________________________________________
[ Info: Pavement time    = 8330.815
[ Info: % time in solver = 8293.543, 100 % of pavement time
[ Info: Single tile time = 20.734
[ Info: _________________________________________________
```