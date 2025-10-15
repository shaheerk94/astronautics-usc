from space_functions import dms_to_decimal

# Dictionary of launch sites
LAUNCH_SITES = {
    "CAPE_CANAVERAL": {
        "name": "Cape Canaveral, Florida",
        "country": "USA",
        "latitude_deg": dms_to_decimal(28, 30, hemisphere='N')
    },
    "VANDENBERG": {
        "name": "Vandenberg Air Force Base, California",
        "country": "USA",
        "latitude_deg": dms_to_decimal(34, 36, hemisphere='N')
    },
    "WALLOPS_ISLAND": {
        "name": "Wallops Island (NASA; Mid-Atl.Reg.SpaceP)",
        "country": "USA",
        "latitude_deg": dms_to_decimal(37, 48, hemisphere='N')
    },
    "KWAJALEIN_ATOLL": {
        "name": "Kwajalein Atoll, RTS (Omelek Island)",
        "country": "USA",
        "latitude_deg": dms_to_decimal(9, 3, hemisphere='N')
    },
    "KODIAK_ISLAND": {
        "name": "Kodiak Island (Alaska Aerospace Dvlp Corp)",
        "country": "USA",
        "latitude_deg": dms_to_decimal(57, 25, hemisphere='N')
    },
    "BAIKONUR": {
        "name": "Baikonur (Tyuratam) Cosmodrome",
        "country": "Russia/Kazakhstan",
        "latitude_deg": dms_to_decimal(45, 54, hemisphere='N')
    },
    "PLESETSK": {
        "name": "Plesetsk Cosmodrome",
        "country": "Russia",
        "latitude_deg": dms_to_decimal(62, 48, hemisphere='N')
    },
    "KAPUSTIN_YAR": {
        "name": "Kapustin Yar (Kap Yar)",
        "country": "Russia",
        "latitude_deg": dms_to_decimal(48, 24, hemisphere='N')
    },
    "VOSTOCHNY": {
        "name": "Vostochny (Svobodny)",
        "country": "Russia",
        "latitude_deg": dms_to_decimal(51, 53, hemisphere='N')
    },
    "KOUROU": {
        "name": "Kourou, European Space Agency (ESA)",
        "country": "French Guiana",
        "latitude_deg": dms_to_decimal(5, 32, hemisphere='N')
    },
    "TANEGASHIMA": {
        "name": "Tanegashima Island",
        "country": "Japan",
        "latitude_deg": dms_to_decimal(30, 14, hemisphere='N')
    },
    "SRIHARIKOTA": {
        "name": "Sriharikota Space Center (Andhra Pradesh)",
        "country": "India",
        "latitude_deg": dms_to_decimal(13, 47, hemisphere='N')
    },
    "JIUQUAN": {
        "name": "Jiuquan Satellite Launch Center (Gobi)",
        "country": "PR China",
        "latitude_deg": dms_to_decimal(40, 25, hemisphere='N')
    },
    "XICHANG": {
        "name": "Xichang Satellite Launch Center (Sichuan)",
        "country": "PR China",
        "latitude_deg": dms_to_decimal(28, 6, hemisphere='N')
    },
    "TAIYUAN": {
        "name": "Taiyuan (Wuzhai) Sat Launch Ctr (Shanxi)",
        "country": "PR China",
        "latitude_deg": dms_to_decimal(38, 40, hemisphere='N')
    },
    "WENCHANG": {
        "name": "Wenchang Sat Launch Center (Hainan)",
        "country": "PR China",
        "latitude_deg": dms_to_decimal(19, 37, hemisphere='N')
    },
    "YAVNE": {
        "name": "Yavne, Palmachim Air Force Base",
        "country": "Israel",
        "latitude_deg": dms_to_decimal(31, 31, hemisphere='N')
    },
    "ALCANTARA": {
        "name": "Alcantara",
        "country": "Brazil",
        "latitude_deg": dms_to_decimal(2, 17, hemisphere='S')
    },
    "MAHIA_PENINSULA": {
        "name": "Mahia Peninsula (Rocket Lab)",
        "country": "New Zealand",
        "latitude_deg": dms_to_decimal(39, 16, hemisphere='S')
    },
    "ARNHEM_LAND": {
        "name": "Arnhem Land (under construction)",
        "country": "Australia",
        "latitude_deg": dms_to_decimal(12, 23, hemisphere='S')
    },
    "WHALERS_WAY": {
        "name": "Whalers Way (Port Lincoln) (under construction)",
        "country": "Australia",
        "latitude_deg": dms_to_decimal(34, 30, hemisphere='S')
    },
    "NINGBO": {
        "name": "Ningbo (Zhejiang) (under construction)",
        "country": "PR China",
        "latitude_deg": dms_to_decimal(29, 50, hemisphere='N')
    }
}
