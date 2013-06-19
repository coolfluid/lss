Coolfluid3 lss plugins
======================

Collection of plugins to solve Linear Systems

**Installation:**

  + Create a plugins directory, and clone the lss sources inside:

```
mkdir -p $PLUGIN_DIR
cd $PLUGIN_DIR
git clone https://github.com/coolfluid/lss.git $PLUGIN_DIR/lss
```

  + Rerun cmake in the coolfluid3 build directory:

```
cd $CF3_BUILD_DIR
cmake .  -DCF3_PLUGIN_DIRS=$PLUGIN_DIR -DCF3_PLUGIN_LSS=ON
```