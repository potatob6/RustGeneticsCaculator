#[macro_export]
// 标准前景色宏定义
macro_rules! black {
    ($n:expr) => {
        format!("\x1b[30m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! red {
    ($n:expr) => {
        format!("\x1b[31m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! green {
    ($n:expr) => {
        format!("\x1b[32m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! yellow {
    ($n:expr) => {
        format!("\x1b[33m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! blue {
    ($n:expr) => {
        format!("\x1b[34m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! magenta {
    ($n:expr) => {
        format!("\x1b[35m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! cyan {
    ($n:expr) => {
        format!("\x1b[36m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! white {
    ($n:expr) => {
        format!("\x1b[37m{}\x1b[0m", $n)
    };
}

#[macro_export]
// 亮色前景色宏定义
macro_rules! bright_black {
    ($n:expr) => {
        format!("\x1b[90m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_red {
    ($n:expr) => {
        format!("\x1b[91m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_green {
    ($n:expr) => {
        format!("\x1b[92m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_yellow {
    ($n:expr) => {
        format!("\x1b[93m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_blue {
    ($n:expr) => {
        format!("\x1b[94m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_magenta {
    ($n:expr) => {
        format!("\x1b[95m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_cyan {
    ($n:expr) => {
        format!("\x1b[96m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_white {
    ($n:expr) => {
        format!("\x1b[97m{}\x1b[0m", $n)
    };
}

// 标准背景色宏定义
#[macro_export]
macro_rules! black_bg {
    ($n:expr) => {
        format!("\x1b[;40m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! red_bg {
    ($n:expr) => {
        format!("\x1b[;41m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! green_bg {
    ($n:expr) => {
        format!("\x1b[;42m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! yellow_bg {
    ($n:expr) => {
        format!("\x1b[;43m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! blue_bg {
    ($n:expr) => {
        format!("\x1b[;44m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! magenta_bg {
    ($n:expr) => {
        format!("\x1b[;45m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! cyan_bg {
    ($n:expr) => {
        format!("\x1b[;46m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! white_bg {
    ($n:expr) => {
        format!("\x1b[;47m{}\x1b[0m", $n)
    };
}

// 明亮背景色宏定义
#[macro_export]
macro_rules! bright_black_bg {
    ($n:expr) => {
        format!("\x1b[;100m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_red_bg {
    ($n:expr) => {
        format!("\x1b[;101m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_green_bg {
    ($n:expr) => {
        format!("\x1b[;102m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_yellow_bg {
    ($n:expr) => {
        format!("\x1b[;103m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_blue_bg {
    ($n:expr) => {
        format!("\x1b[;104m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_magenta_bg {
    ($n:expr) => {
        format!("\x1b[;105m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_cyan_bg {
    ($n:expr) => {
        format!("\x1b[;106m{}\x1b[0m", $n)
    };
}

#[macro_export]
macro_rules! bright_white_bg {
    ($n:expr) => {
        format!("\x1b[;107m{}\x1b[0m", $n)
    };
}