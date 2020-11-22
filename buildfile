cxx.std = latest
using cxx

hxx{*}: extension = hpp
cxx{*}: extension = cpp

import libs = glfw3%lib{glfw3}
import libs += glbinding%lib{glbinding}
import libs += glm%lib{glm}

exe{ray}: {hxx cxx}{**} $libs