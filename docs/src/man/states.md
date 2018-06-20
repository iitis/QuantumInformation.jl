```@setup QI
Pkg.add("QI")
importall QI
```

# States and channels

## Generalized transpositions
## States
```@repl QI
ket(0,2)  
ket(1,8)
ket(1,8,sparse=true)
```

```@repl QI
bra(0,2)  
bra(1,4)
bra(1,4,sparse=true)
```

```@repl QI
ketbra(0,1,2)  
ketbra(1,3,4)
ketbra(1,3,4,sparse=true)
```

```@repl QI
proj(ket(1,2))
proj(ket(3,4,sparse=true))
```

```@repl QI
res(ketbra(0,1,2))
```

```@repl QI
res(ketbra(0,1,2))
```

```@repl QI
unres(res(ketbra(0,1,2)))
```

## Channels
