use rug::{Float,ops::Pow};
// use rug::{ops::Pow, Float};
fn rat(prec: u32, num: &str, den: &str) -> Float{
    let num = Float::with_val(prec, Float::parse(num).unwrap());
    let den = Float::with_val(prec, Float::parse(den).unwrap());

    // println!("n/d = {}", num/den);
    num/den
}


// Generators initialization
pub struct Generators{
    lambda1 : f64,
    lambda2 : f64,
    lambda3 : f64,
    lambda4 : f64,
    lambda5 : f64,
    delta   : f64,
}

impl Generators {
    pub fn new(prec: u32) -> Generators {
        Generators {

            // TODO increase precision
            lambda1 : rat(prec, "18828",
                          "20605").sqrt().to_f64(),
            lambda2 : rat(prec, "51249816",
                          "310826425").sqrt().to_f64(),
            lambda3 : rat(prec, "16769730234298044",
                          "20923428441426875").sqrt().to_f64(),
            lambda4 : rat(prec, "7291912559866455728218951411427664",
                          "10152569850187822006207489090234375").sqrt().to_f64(),
            lambda5 : rat(prec, "4650226699147357157251468602113895310677921530735338408908488478132",
                          "6022803457906407975449078147096424232221963052245142351118896484375").sqrt().to_f64(),
            delta   : rat(prec, "4707","10000").sqrt().to_f64(),
        }
    }

    pub fn get_l1(&self) -> f64 { self.lambda1 }
    pub fn get_l2(&self) -> f64 { self.lambda2 }
    pub fn get_l3(&self) -> f64 { self.lambda3 }
    pub fn get_l4(&self) -> f64 { self.lambda4 }
    pub fn get_l5(&self) -> f64 { self.lambda5 }
    pub fn get_de(&self) -> f64 { self.delta   }

    // Squares for difference
    pub fn get_l1s(&self) -> f64 { self.lambda1*self.lambda1 }
    pub fn get_l2s(&self) -> f64 { self.lambda2*self.lambda2 }
    pub fn get_l3s(&self) -> f64 { self.lambda3*self.lambda3 }
    pub fn get_l4s(&self) -> f64 { self.lambda4*self.lambda4 }
    pub fn get_l5s(&self) -> f64 { self.lambda5*self.lambda5 }

    pub fn info (&self) {
        println!("lambda[1] = {}", self.lambda1);
        println!("lambda[2] = {}", self.lambda2);
        println!("lambda[3] = {}", self.lambda3);
        println!("lambda[4] = {}", self.lambda4);
        println!("lambda[5] = {}", self.lambda5);
        println!("delta     = {}", self.delta);
    }
}

// Rule of the order 11
pub struct Rule11weights{
    n_dim  : u32,
    w00000: f64,
    w10000: f64,
    w11000: f64,
    w11100: f64,
    w11110: f64,
    w20000: f64,
    w21000: f64,
    w21100: f64,
    w22000: f64,
    w30000: f64,
    w31000: f64,
    w40000: f64,
    w_delta: f64,
}




impl Rule11weights {
    pub fn new(prec: u32, n: usize) -> Rule11weights {
        let n_dim:u32 = n as u32;
        Rule11weights {
            n_dim,

            w00000: ((rat(prec, "461470916330619094601", "561470916330619094601") -
                     rat(prec, "178776407247243768611164181152393408208408973929188273376004254015863055",
                         "206645117428890102583802582449510927672430049005571921797238974120104704")*n_dim +
                     rat(prec, "292815442456570754577795780296723491558075",
                         "1067131356556584241707053820923415207090816")*n_dim*n_dim -
                     rat(prec, "11036015794961627400432125", "586876352352746947393506048")*n_dim*n_dim*n_dim + 
                      rat(prec, "742836251796902875625", "3449677309935323717228544")*n_dim*n_dim*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w10000: ((rat(prec, "10772827971254467627487919796646539771484336444230118029252981549368377815",
                         "69526866183422813820487555648490525652524987918091268761829655499440313024") -
                     rat(prec, "7823019915955833629860880942551959935934415075",
                         "71219833948998491067359722476320807399165058816")*n_dim +
                     rat(prec, "21018950990011229581625", "1018882556167963450335948")*n_dim*n_dim  - 
                      rat(prec, "742836251796902875625", "1724838654967661858614272")*n_dim*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w11000: ((rat(prec, "50606762189236105532484922234975","2961524336600664120036854144215296") -
                     rat(prec, "5106538565366996918255875", "391250901568497964929004032")*n_dim + 
                      rat(prec, "742836251796902875625", "1149892436645107905742848")*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(),


            w11100: ((rat(prec, "786847887791947713875", "287473109161276976435712") -
                      rat(prec, "742836251796902875625", "1149892436645107905742848")*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w11110: (rat(prec, "742836251796902875625", "2299784873290215811485696")*u32::pow(2,n_dim)).to_f64(),

            w20000: ((rat(prec, "789632406937983589159627892938968218048307729053605005635326529375",
                         "2272983234278990247670015265977739148998772364018378608399912647424") -
                     rat(prec, "15736021731215429025618229864025553125", "91412413588103446613731884361225330176")*n_dim +
                      rat(prec, "4289984865961615259359375", "521667868757997286572005376")*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(), 

            w21000: ((rat(prec, "847193900018576642279449061303125", "22388541167794133385680108832041472") -
                      rat(prec, "4289984865961615259359375", "521667868757997286572005376")*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w21100: (rat(prec, "4289984865961615259359375", "1043335737515994573144010752")*u32::pow(2,n_dim)).to_f64(),

            w22000: (rat(prec, "24775271946759384235234765625", "473326646253089538016332877824")*u32::pow(2,n_dim)).to_f64(),

            w30000: ((rat(prec, "-70934649605204746122516339068114323575860845274075071699750436295559326171875",
                         "1659483269533274798949736283963096149092084677977802804132444108445406644778752") -
                     rat(prec, "1037385528652157530660005994972455364403076171875",
                         "49290927916765393329088241190712613882690601059328")*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w31000: (rat(prec, "1037385528652157530660005994972455364403076171875",
                         "98581855833530786658176482381425227765381202118656")*u32::pow(2,n_dim)).to_f64(),

            w40000: (rat(prec, "12310921778223316311723101101343809674162717353025847744054643753256\
                               4812928669370188649233161387562544199940585531294345855712890625",
                        "10286772307809702183778730306541235862931879876979919480490032496827154255\
                         55748479498967391959526098968217357319751238154101146595328")*u32::pow(2,n_dim)).to_f64(),

            w_delta: (rat(prec, "10000", "14121").to_f64()).powi(5), // FIXME
        }
    }


    pub fn get_w00000(&self) -> f64 { self.w00000 }
    pub fn get_w10000(&self) -> f64 { self.w10000 }
    pub fn get_w11000(&self) -> f64 { self.w11000 }
    pub fn get_w11100(&self) -> f64 { self.w11100 }
    pub fn get_w11110(&self) -> f64 { self.w11110 }
    pub fn get_w20000(&self) -> f64 { self.w20000 }
    pub fn get_w21000(&self) -> f64 { self.w21000 }
    pub fn get_w21100(&self) -> f64 { self.w21100 }
    pub fn get_w22000(&self) -> f64 { self.w22000 }
    pub fn get_w30000(&self) -> f64 { self.w30000 }
    pub fn get_w31000(&self) -> f64 { self.w31000 }
    pub fn get_w40000(&self) -> f64 { self.w40000 }
    pub fn get_w_delta(&self) -> f64 { self.w_delta }

    pub fn info (&self) {
        println!("d      = {}", self.n_dim);
        println!("w00000 = {}", self.w00000);
        println!("w00000 = {}", self.w10000);
        println!("w00000 = {}", self.w11000);
        println!("w00000 = {}", self.w11100);
        println!("w00000 = {}", self.w11110);
        println!("w00000 = {}", self.w20000);
        println!("w00000 = {}", self.w21000);
        println!("w00000 = {}", self.w21100);
        println!("w00000 = {}", self.w22000);
        println!("w00000 = {}", self.w30000);
        println!("w00000 = {}", self.w31000);
        println!("w40000 = {}", self.w40000);
        println!("w_delta = {}", self.w_delta);
    }
}




// Rule of the order 13
pub struct Rule13weights{
    n_dim  : u32,
    w00000: f64,
    w10000: f64,
    w11000: f64,
    w11100: f64,
    w11110: f64,
    w11111: f64,
    w20000: f64,
    w21000: f64,
    w21100: f64,
    w21110: f64,
    w22000: f64,
    w22100: f64,
    w30000: f64,
    w31000: f64,
    w31100: f64,
    w32000: f64,
    w40000: f64,
    w41000: f64,
    w50000: f64,
    w_delta: f64,
}




impl Rule13weights {
    pub fn new(prec: u32, n: usize) -> Rule13weights {
        let n_dim:u32 = n as u32;
        Rule13weights {
            n_dim,

            w00000: ((rat(prec, "6928530809504672234860721", "7928530809504672234860721") -
                     rat(prec, "743347207984408938393222010616297422021070021340478815\
                                413981089859624110276163422146328513037894196223861577\
                                408130755169407116044983037445",
                         "720709981737199219527198605885667582383725953904867538164369\
                          058784783393444148884346667637746826334458637834457064200782921566241090749696")*n_dim +
                     rat(prec, "20567895888045881470775951133655472423774547011880265671985900726786297025",
                         "48033509517915343845035000276041873410080406946628495599980437095473226752")*n_dim*n_dim -
                     rat(prec, "1893915787204371571584772357973364273222384125",
                         "30137923771871052154290614010519092278658825472")*n_dim*n_dim*n_dim +
                     rat(prec, "47823097396099162653325785625", "22099415924195039051049863743488")*n_dim*n_dim*n_dim*n_dim - 
                      rat(prec, "3061228193655036750450625", "194851573174386824843937079296")*n_dim*n_dim*n_dim*n_dim*n_dim )*u32::pow(2,n_dim)).to_f64(),

            w10000: ((rat(prec, "114605223550967205439295880384886047911791381783068794694\
                                698704561189995016719858039882687354523517108258646726933\
                                077607277254723297986455815789124983990216999698215458085",
                         "868515400220762225023735912224912599135810857885017247225104436\
                          353079257645573141585284972149138773631862661459963159126058400\
                          907940361820962880367599477421357990956911616") -
                     rat(prec,"4675070889152040578633823057748352015904886876592465975833\
                               07136543946814099636412090090071228952868719730475",
                         "301385339413572177540504609567886178533260781544545125254935265\
                          2171139258163028238943966747970095029032531968")*n_dim +
                     rat(prec, "530495365893417907651636661867970203289037681805777875",
                         "10950010156310178154439488148167517208535943453847533568")*n_dim*n_dim -
                     rat(prec, "225850145218696295377333898125", "66298247772585117153149591230464")*n_dim*n_dim*n_dim +
                      rat(prec, "15306140968275183752253125", "389703146348773649687874158592")*n_dim*n_dim*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w11000: ((rat(prec, "10104799241174402564154880744835326314756209136749938107631076703184357075",
                         "317653927809144561662737126364906483131701158097991363321458081471357004032") -
                     rat(prec, "18951195107154009644032269370134412619242410997625",
                         "1005695275193807692362186641088126121283609795540736")*n_dim +
                     rat(prec, "41190426515199403708678270625", "11049707962097519525524931871744")*n_dim*n_dim -
                      rat(prec, "15306140968275183752253125", "194851573174386824843937079296")*n_dim*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(), 

            w11100: ((rat(prec, "5473008354217394095939886767310347375", "2007344887542622945873940033782281351168") -
                     rat(prec, "11519251878099881588010251875", "4910981316487786455788858609664")*n_dim +
                      rat(prec, "15306140968275183752253125", "129901048782924549895958052864")*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(), 

            w11110: ((rat(prec, "15644683227406412104875625", "32475262195731137473989513216") -
                      rat(prec, "15306140968275183752253125",  "129901048782924549895958052864")*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w11111: (rat(prec, "15306140968275183752253125", "259802097565849099791916105728")*u32::pow(2,n_dim)).to_f64(),

            w20000: ((rat(prec, "34964101366258726399629871621474546415855151193087537616126821605371069876262\
                                7959593296943640100071043797485831844737926012175472831357017729375",
                         "79328045908373383416602035381901663252328673978322535726088259201596399371\
                          4764097613001826439291972497098282534227450893164406174623779612408064") -
                     rat(prec, "37430231584033668780776861836117735439817251687186356935586314648552636987334375", 
                         "124868663307498822856752264210518703811248089271359486164513174988053120810753536")*n_dim + 
                     rat(prec, "20888800779449717241396562047214108796875",  "430278230759202923210835979688287629138432")*n_dim*n_dim + 
                      rat(prec, "-88395138163139082419099921875", "88397663696780156204199454973952")*n_dim*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w21000: ((rat(prec, "329580969780830729612470764445553569266275233806752006433909924074882678125", 
                         "5659658910238965373515875032740017602054440569344928204579657367782266048512") - 
                     rat(prec, "319639700845602012538550142002852261828125", "10326677538220870157060063512518903099322368")*n_dim +
                      rat(prec, "88395138163139082419099921875", "58931775797853437469466303315968")*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(), 

            w21100: ((rat(prec, "50678979788541811639817338307335015625", "7587566155930102980940531603614183026688") -
                      rat(prec, "88395138163139082419099921875", "58931775797853437469466303315968")*n_dim)*u32::pow(2,n_dim)).to_f64(), 

            w21110: (rat(prec, "88395138163139082419099921875", "117863551595706874938932606631936")*u32::pow(2,n_dim)).to_f64(), 

            w22000: ((rat(prec, "105544490536207355000531777199548454296875", "1147408615357874461895562612502100344369152") -
                      rat(prec, "510494478462977112167012345703125", "26735382286959509465314546271010816")*n_dim)*u32::pow(2,n_dim)).to_f64(), 

            w22100: (rat(prec, "510494478462977112167012345703125", "53470764573919018930629092542021632")*u32::pow(2,n_dim)).to_f64(), 

            w30000: ((rat(prec, "303303336986482428831411139090668345284780166846006526378126448204989966794\
                                5338129639962460610772733791808491484550108334904587137377786602294921875",
                         "3596575990445299293965504129417615247070179376089292142943738322880952738058\
                          088530052588312355225816582396620141932718309845524429720507407952267264") -
                     rat(prec, "583943925500194570596058038432163509635078008660102711302768170928439722900390625", 
                         "20637760790744658899201887357833639512564403494972493520036835765132525258453870592")*n_dim +
                     rat(prec, "21375328817877705919249423526407442783525384521484375",
                         "5568297544901152953600440430832422565099791820470165504")*n_dim*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w31000: ((rat(prec, "-380609151806305191962294716142902945105148768315681002510192291376075439453125", 
                         "45491023050869931445705850164218162040920801237999618339537477807051855823190016") -
                     rat(prec, "21375328817877705919249423526407442783525384521484375",
                         "5568297544901152953600440430832422565099791820470165504")*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w31100: (rat(prec, "21375328817877705919249423526407442783525384521484375",
                         "11136595089802305907200880861664845130199583640940331008")*u32::pow(2,n_dim)).to_f64(),

            w32000: (rat(prec, "123445560056915451191546345478520420265625915679931640625",
                         "5052301972406979446566799617575284740733877778439930167296")*u32::pow(2,n_dim)).to_f64(),

            w40000: ((rat(prec, "1773115147467478873651784244790993435647570577211102447704944703\
                                186191528218882929869570040379807543393674873771910873139389908\
                                8508644541623021959340432904442977817831466439715586602687835693359375",
                         "278128092409622915464524677580623060271292281178977037749589195847539\
                          123841920361805828456960981295895111114317773383685984381737828082283\
                          75479726062573728793036655972628153015788984323690863300608") -
                     rat(prec, "6155460889111658155861550550671904837081358676512923872027321\
                                87662824064643346850943246165806937812720999702927656471729278564453125",
                         "14099443024370861881789803509698499550639269666860804948701747040785\
                          852486554452102892422753565123070643239313479453757800594313101312")*n_dim)*u32::pow(2,n_dim)).to_f64(),

            w41000: (rat(prec, "615546088911165815586155055067190483708135867651292387202732187\
                               662824064643346850943246165806937812720999702927656471729278564453125",
                        "281988860487417237635796070193969991012785393337216098974034940815717\
                         04973108904205784845507130246141286478626958907515601188626202624")*u32::pow(2,n_dim)).to_f64(),

            w50000: (rat(prec, "-7665293024343553319250239729049240774493191908857924342705617\
                               86673657822213862992959630370078613151039226836021436756503717\
                               43623708053039473017541376891935290029388776647794020498857987\
                               11031049271820895918931052152724231796508155431278315041644951\
                               95541510633090644720486123488566000104148245863910915431915782392024993896484375",
                        "59364615520007290964320282556730475794503997684545374969991463593755\
                         10906949371731021249665468840756829806883780853165380992852995196238\
                         60209388275765879788321248877982461591141134025621223715632640756863\
                         34156814373716395262095137572920740443778798303205490353342110809838\
                         2247366208522693568318641632774369089945244275844366336")*u32::pow(2,n_dim)).to_f64(),

            w_delta: (rat(prec, "10000", "14121").to_f64()).powi(6),
        }
    }

    pub fn get_w00000(&self) -> f64 { self.w00000 }
    pub fn get_w10000(&self) -> f64 { self.w10000 }
    pub fn get_w11000(&self) -> f64 { self.w11000 }
    pub fn get_w11100(&self) -> f64 { self.w11100 }
    pub fn get_w11110(&self) -> f64 { self.w11110 }
    pub fn get_w11111(&self) -> f64 { self.w11111 }
    pub fn get_w20000(&self) -> f64 { self.w20000 }
    pub fn get_w21000(&self) -> f64 { self.w21000 }
    pub fn get_w21100(&self) -> f64 { self.w21100 }
    pub fn get_w21110(&self) -> f64 { self.w21110 }
    pub fn get_w22000(&self) -> f64 { self.w22000 }
    pub fn get_w22100(&self) -> f64 { self.w22100 }
    pub fn get_w30000(&self) -> f64 { self.w30000 }
    pub fn get_w31000(&self) -> f64 { self.w31000 }
    pub fn get_w31100(&self) -> f64 { self.w31100 }
    pub fn get_w32000(&self) -> f64 { self.w32000 }
    pub fn get_w40000(&self) -> f64 { self.w40000 }
    pub fn get_w41000(&self) -> f64 { self.w41000 }
    pub fn get_w50000(&self) -> f64 { self.w50000 }
    pub fn get_w_delta(&self) -> f64 { self.w_delta }


    pub fn info (&self) {
        println!("d      = {}", self.n_dim);
        println!("w00000 = {}", self.w00000);
        println!("w10000 = {}", self.w10000);
        println!("w11000 = {}", self.w11000);
        println!("w11100 = {}", self.w11100);
        println!("w11110 = {}", self.w11110);
        println!("w11111 = {}", self.w11111);
        println!("w20000 = {}", self.w20000);
        println!("w21000 = {}", self.w21000);
        println!("w21100 = {}", self.w21100);
        println!("w21110 = {}", self.w21110);
        println!("w22000 = {}", self.w22000);
        println!("w22100 = {}", self.w22100);
        println!("w30000 = {}", self.w30000);
        println!("w31000 = {}", self.w31000);
        println!("w31100 = {}", self.w31100);
        println!("w32000 = {}", self.w32000);
        println!("w40000 = {}", self.w40000);
        println!("w41000 = {}", self.w41000);
        println!("w50000 = {}", self.w50000);
        println!("w_delta = {}", self.w_delta);
    }
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn rule11_weights() {
        assert_eq!(2 + 2, 4);
        
        let w11d5 = Rule11weights::new(100,5);
        w11d5.info();
        
        // w11.w00000 = Float::with_val(5);
        println!("w000 = {}", w11d5.w00000)
    }
    
}
