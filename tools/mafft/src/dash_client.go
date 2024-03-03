package main

import(
    "bufio"
    "bytes"
    "encoding/json"
    "flag"
    "fmt"
    "io"
    "io/ioutil"
    "net/http"
    "os"
    "strings"
    "strconv"
    "sort"
    "time"
)

var VersionNumber = 1.1
var Debug = false
var NumberOfQueryTries = 12
var HTTPClient = &http.Client{}
var SequenceLimit = 3000
var AlignmentLimit = 10000

//URL's
func DASHDomainAlignmentURL(dash_url string) string {
    return dash_url + "domain_alignments?format=JSON"
}
func DASHChainDomainsURL(dash_url string) string {
    return dash_url + "domains?format=JSON&filter=pdbid=%s"
}
func DASHChainURL(dash_url string) string {
    return dash_url + "chains?format=JSON&filter=pdbid=%s"
}
func DASHChainSearchURL(dash_url string) string {
    return dash_url + "chain_search_sequence?limit=5"
}

//Data types
type Sequence struct {
    Label string
    Sequence string
}
type DASHInput struct {
    PDBID string
    FullID string
    Hat3Index int
    Start int
    End int
    Sequence string
    Domains []RESTDomain
}
type RESTChain struct {
    StatusCode int
    StatusMessage string
    PDBID string
    DepositionDate string
    Sequence string
    Length int
}
type RESTDomain struct {
    StatusCode int
    StatusMessage string
    DomainID string
    PDBID string
    Length int
    Sequence string
    Segments string
    ResidueNumbers string
    Start int
    End int
    SliceStart int
    SliceEnd int
    ResidueNumberInts []int
}
type RESTAlignment struct {
    StatusCode int
    StatusMessage string
    SCORE int
    ID1 string
    ID2 string
    PRIMS1 string
    SECOS1 string
    PRIMS2 string
    SECOS2 string
    EQUIVALENCE string
    LOWSIMILARITY bool
}
type RESTSearch struct {
    StatusCode int
    StatusMessage string
    ID string
    Start int
    End int
}
func(alignment *RESTAlignment) Reverse() {
    new_alignment := *alignment
    new_alignment.ID1 = alignment.ID2
    new_alignment.ID2 = alignment.ID1
    new_alignment.PRIMS1 = alignment.PRIMS2
    new_alignment.PRIMS2 = alignment.PRIMS1
    new_alignment.SECOS1 = alignment.SECOS2
    new_alignment.SECOS2 = alignment.SECOS1
    *alignment = new_alignment
}

//Utility Functions
func fatal(object interface{}) {
    if Debug {
        panic(fmt.Sprint(object))
    } else {
        fmt.Fprintln(os.Stderr, "----")
        non_fatal(object)
        os.Exit(1)
    }
}

func non_fatal(object interface{}) {
    fmt.Fprintln(os.Stderr, fmt.Sprint(object))
}

func check(err error) {
	if err != nil {
        fatal(err.Error())
	}
}

func http_query(method string, url string, body io.Reader) *http.Response {
    //A simple request will follow redirects by default!
    retry_count := NumberOfQueryTries
    var response *http.Response
    status_code := 404
    for status_code != 200 {
        if retry_count != NumberOfQueryTries {
            fmt.Println("Retrying DASH request...")
            time.Sleep(10*time.Second)
        }
        if retry_count == 0 {
            break
        }
        request, err := http.NewRequest(method, url, body)
        check(err)
        response, err = HTTPClient.Do(request)
        if err != nil {
            status_code = 404
            retry_count += -1
            continue
        }
        status_code = response.StatusCode
    }

    if status_code != 200 {
        error_message := fmt.Sprintf(
            "Client was unable to connect to DASH server after %d retries.", NumberOfQueryTries)
        error_message +=
            "\nPlease check https://sysimm.org for information about possible maintenance."
        error_message +=
            "\nIf there is no scheduled maintenance occuring right now you may submit a bug report by contacting us at https://sysimm.org"
        fatal(error_message)
    }

    return response
}

func parse_residue_numbers(residue_numbers_string string) []int {
    starts_and_ends := strings.Split(residue_numbers_string, "; ")
    residue_numbers := make([]int, 0, 3000)
    for _, start_and_end := range(starts_and_ends) {
        start_and_end_split := strings.Split(start_and_end, "-")
        start, err := strconv.Atoi(start_and_end_split[0])
        check(err)
        end, err := strconv.Atoi(start_and_end_split[1])
        check(err)
        for residue_number := start; residue_number <= end; residue_number++ {
            residue_numbers = append(residue_numbers, residue_number)
        }
    }
    return residue_numbers
}

//BLOSUM62 Matrix
var BLOSUM62Max = 11.0
type BLOSUMRow map[byte]float64
type BLOSUMMatrix map[byte]BLOSUMRow
var BLOSUM62 = BLOSUMMatrix{
    'A':BLOSUMRow{
        'A':4, 'R':-1, 'N':-2, 'D':-2, 'C':0, 'Q':-1, 'E':-1, 'G':0, 'H':-2, 'I':-1,
        'L':-1, 'K':-1, 'M':-1, 'F':-2, 'P':-1, 'S':1, 'T':0, 'W':-3, 'Y':-2, 'V':0,
        'B':-2, 'Z':-1, },
    'R':BLOSUMRow{
        'A':-1, 'R':5, 'N':0, 'D':-2, 'C':-3, 'Q':1, 'E':0, 'G':-2, 'H':0, 'I':-3,
        'L':-2, 'K':2, 'M':-1, 'F':-3, 'P':-2, 'S':-1, 'T':-1, 'W':-3, 'Y':-2, 'V':-3,
        'B':-1, 'Z':0, },
    'N':BLOSUMRow{
        'A':-2, 'R':0, 'N':6, 'D':1, 'C':-3, 'Q':0, 'E':0, 'G':0, 'H':1, 'I':-3,
        'L':-3, 'K':0, 'M':-2, 'F':-3, 'P':-2, 'S':1, 'T':0, 'W':-4, 'Y':-2, 'V':-3,
        'B':3, 'Z':0, },
    'D':BLOSUMRow{
        'A':-2, 'R':-2, 'N':1, 'D':6, 'C':-3, 'Q':0, 'E':2, 'G':-1, 'H':-1, 'I':-3,
        'L':-4, 'K':-1, 'M':-3, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3,
        'B':4, 'Z':1, },
    'C':BLOSUMRow{
        'A':0, 'R':-3, 'N':-3, 'D':-3, 'C':9, 'Q':-3, 'E':-4, 'G':-3, 'H':-3, 'I':-1,
        'L':-1, 'K':-3, 'M':-1, 'F':-2, 'P':-3, 'S':-1, 'T':-1, 'W':-2, 'Y':-2, 'V':-1,
        'B':-3, 'Z':-3, },
    'Q':BLOSUMRow{
        'A':-1, 'R':1, 'N':0, 'D':0, 'C':-3, 'Q':5, 'E':2, 'G':-2, 'H':0, 'I':-3,
        'L':-2, 'K':1, 'M':0, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-2, 'Y':-1, 'V':-2,
        'B':0, 'Z':3, },
    'E':BLOSUMRow{
        'A':-1, 'R':0, 'N':0, 'D':2, 'C':-4, 'Q':2, 'E':5, 'G':-2, 'H':0, 'I':-3,
        'L':-3, 'K':1, 'M':-2, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2,
        'B':1, 'Z':4, },
    'G':BLOSUMRow{
        'A':0, 'R':-2, 'N':0, 'D':-1, 'C':-3, 'Q':-2, 'E':-2, 'G':6, 'H':-2, 'I':-4,
        'L':-4, 'K':-2, 'M':-3, 'F':-3, 'P':-2, 'S':0, 'T':-2, 'W':-2, 'Y':-3, 'V':-3,
        'B':-1, 'Z':-2, },
    'H':BLOSUMRow{
        'A':-2, 'R':0, 'N':1, 'D':-1, 'C':-3, 'Q':0, 'E':0, 'G':-2, 'H':8, 'I':-3,
        'L':-3, 'K':-1, 'M':-2, 'F':-1, 'P':-2, 'S':-1, 'T':-2, 'W':-2, 'Y':2, 'V':-3,
        'B':0, 'Z':0, },
    'I':BLOSUMRow{
        'A':-1, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-3, 'E':-3, 'G':-4, 'H':-3, 'I':4,
        'L':2, 'K':-3, 'M':1, 'F':0, 'P':-3, 'S':-2, 'T':-1, 'W':-3, 'Y':-1, 'V':3,
        'B':-3, 'Z':-3, },
    'L':BLOSUMRow{
        'A':-1, 'R':-2, 'N':-3, 'D':-4, 'C':-1, 'Q':-2, 'E':-3, 'G':-4, 'H':-3, 'I':2,
        'L':4, 'K':-2, 'M':2, 'F':0, 'P':-3, 'S':-2, 'T':-1, 'W':-2, 'Y':-1, 'V':1,
        'B':-4, 'Z':-3, },
    'K':BLOSUMRow{
        'A':-1, 'R':2, 'N':0, 'D':-1, 'C':-3, 'Q':1, 'E':1, 'G':-2, 'H':-1, 'I':-3,
        'L':-2, 'K':5, 'M':-1, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2,
        'B':0, 'Z':1, },
    'M':BLOSUMRow{
        'A':-1, 'R':-1, 'N':-2, 'D':-3, 'C':-1, 'Q':0, 'E':-2, 'G':-3, 'H':-2, 'I':1,
        'L':2, 'K':-1, 'M':5, 'F':0, 'P':-2, 'S':-1, 'T':-1, 'W':-1, 'Y':-1, 'V':1,
        'B':-3, 'Z':-1, },
    'F':BLOSUMRow{
        'A':-2, 'R':-3, 'N':-3, 'D':-3, 'C':-2, 'Q':-3, 'E':-3, 'G':-3, 'H':-1, 'I':0,
        'L':0, 'K':-3, 'M':0, 'F':6, 'P':-4, 'S':-2, 'T':-2, 'W':1, 'Y':3, 'V':-1,
        'B':-3, 'Z':-3, },
    'P':BLOSUMRow{
        'A':-1, 'R':-2, 'N':-2, 'D':-1, 'C':-3, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-3,
        'L':-3, 'K':-1, 'M':-2, 'F':-4, 'P':7, 'S':-1, 'T':-1, 'W':-4, 'Y':-3, 'V':-2,
        'B':-2, 'Z':-1, },
    'S':BLOSUMRow{
        'A':1, 'R':-1, 'N':1, 'D':0, 'C':-1, 'Q':0, 'E':0, 'G':0, 'H':-1, 'I':-2,
        'L':-2, 'K':0, 'M':-1, 'F':-2, 'P':-1, 'S':4, 'T':1, 'W':-3, 'Y':-2, 'V':-2,
        'B':0, 'Z':0, },
    'T':BLOSUMRow{
        'A':0, 'R':-1, 'N':0, 'D':-1, 'C':-1, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-1,
        'L':-1, 'K':-1, 'M':-1, 'F':-2, 'P':-1, 'S':1, 'T':5, 'W':-2, 'Y':-2, 'V':0,
        'B':-1, 'Z':-1, },
    'W':BLOSUMRow{
        'A':-3, 'R':-3, 'N':-4, 'D':-4, 'C':-2, 'Q':-2, 'E':-3, 'G':-2, 'H':-2, 'I':-3,
        'L':-2, 'K':-3, 'M':-1, 'F':1, 'P':-4, 'S':-3, 'T':-2, 'W':11, 'Y':2, 'V':-3,
        'B':-4, 'Z':-3, },
    'Y':BLOSUMRow{
        'A':-2, 'R':-2, 'N':-2, 'D':-3, 'C':-2, 'Q':-1, 'E':-2, 'G':-3, 'H':2, 'I':-1,
        'L':-1, 'K':-2, 'M':-1, 'F':3, 'P':-3, 'S':-2, 'T':-2, 'W':2, 'Y':7, 'V':-1,
        'B':-3, 'Z':-2, },
    'V':BLOSUMRow{
        'A':0, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-2, 'E':-2, 'G':-3, 'H':-3, 'I':3,
        'L':1, 'K':-2, 'M':1, 'F':-1, 'P':-2, 'S':-2, 'T':0, 'W':-3, 'Y':-1, 'V':4,
        'B':-3, 'Z':-2, },
    'B':BLOSUMRow{
        'A':-2, 'R':-1, 'N':3, 'D':4, 'C':-3, 'Q':0, 'E':1, 'G':-1, 'H':0, 'I':-3,
        'L':-4, 'K':0, 'M':-3, 'F':-3, 'P':-2, 'S':0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3,
        'B':4, 'Z':1, },
    'Z':BLOSUMRow{
        'A':-1, 'R':0, 'N':0, 'D':1, 'C':-3, 'Q':3, 'E':4, 'G':-2, 'H':0, 'I':-3,
        'L':-3, 'K':1, 'M':-1, 'F':-3, 'P':-1, 'S':0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2,
        'B':1, 'Z':4, },
}

//Realign
type FloatRow []float64
type FloatMatrix []FloatRow
type IntRow []int
type IntMatrix []IntRow
func InitializeFloatMatrix(size_a int, size_b int) *FloatMatrix {
    float_matrix := make(FloatMatrix, size_a)
    for i, _ := range(float_matrix) {
        float_matrix[i] = make(FloatRow, size_b)
    }
    return &float_matrix
}
func InitializeIntMatrix(size_a int, size_b int) *IntMatrix {
    int_matrix := make(IntMatrix, size_a)
    for i, _ := range(int_matrix) {
        int_matrix[i] = make(IntRow, size_b)
    }
    return &int_matrix
}
func AlignMatrix(equivalence_matrix *FloatMatrix) (FloatMatrix, int, float64) {
    number_of_residues_a := len(*equivalence_matrix)
    number_of_residues_b := len((*equivalence_matrix)[0])
    //Constants
    bog := 0.0
    beg := 0.0
    iog := 5.0
    ieg := 1.0

    //Initialize matrices
    D := InitializeFloatMatrix(number_of_residues_a+1, number_of_residues_b+1)
    E := InitializeFloatMatrix(number_of_residues_a+1, number_of_residues_b+1)
    F := InitializeFloatMatrix(number_of_residues_a+1, number_of_residues_b+1)
    P1 := InitializeIntMatrix(number_of_residues_a+1, number_of_residues_b+1)
    P2 := InitializeIntMatrix(number_of_residues_a+1, number_of_residues_b+1)
    pE := InitializeIntMatrix(number_of_residues_a+1, number_of_residues_b+1)
    pF := InitializeIntMatrix(number_of_residues_a+1, number_of_residues_b+1)

    //Fill matrices
    for i := 0; i <= number_of_residues_a; i++ {
        if i != number_of_residues_a {
            (*D)[i][number_of_residues_b] =
                -1 * (bog + float64(number_of_residues_a - 1 - i) * beg)
            (*P1)[i][number_of_residues_b] = number_of_residues_a
            (*P2)[i][number_of_residues_b] = number_of_residues_a
        }
        (*E)[i][number_of_residues_b] = -2000.0
        (*F)[i][number_of_residues_b] = -2000.0
        (*pE)[i][number_of_residues_b] = number_of_residues_a
        (*pF)[i][number_of_residues_b] = number_of_residues_a
    }

    for i := 0; i <= number_of_residues_b; i++ {
        if i != number_of_residues_b {
            (*D)[number_of_residues_a][i] =
                -1 * (bog + float64(number_of_residues_b - 1 - i) * beg)
            (*P1)[number_of_residues_a][i] = number_of_residues_b
            (*P2)[number_of_residues_a][i] = number_of_residues_b
        }
        (*E)[number_of_residues_a][i] = -2000.0
        (*F)[number_of_residues_a][i] = -2000.0
        (*pE)[number_of_residues_a][i] = number_of_residues_b
        (*pF)[number_of_residues_a][i] = number_of_residues_b
    }

    for i := 0; i < number_of_residues_a; i++ {
        for j := 0; j < number_of_residues_b; j++ {
            similarity := (*equivalence_matrix)[i][j]
            (*D)[i][j] = similarity
        }
    }

    //Solv_Rec
    for i := number_of_residues_a - 1; i >= 0; i-- {
        gp1o := iog
        gp1e := ieg
        if i == 0 {
            gp1o = bog
            gp1e = beg
        }
        for j := number_of_residues_b - 1; j >= 0; j-- {
            gp2o := iog
            gp2e := ieg
            if j == 0 {
                gp2o = bog
                gp2e = beg
            }

            //Determine E
            d1 := (*E)[i+1][j] - gp2e
            d2 := (*D)[i+1][j] - gp2o
            if d1 > d2 {
                (*E)[i][j] = d1
                (*pE)[i][j] = (*pE)[i+1][j]
            } else {
                (*E)[i][j] = d2
                (*pE)[i][j] = i+1
            }
            //Determine F
            d1 = (*F)[i][j+1] - gp1e
            d2 = (*D)[i][j+1] - gp1o
            if d1 > d2 {
                (*F)[i][j] = d1;
                (*pF)[i][j] = (*pF)[i][j+1]
            } else {
                (*F)[i][j] = d2;
                (*pF)[i][j] = j+1
            }
            //Determine D
            Mx := 0.0
            if (*E)[i][j] > (*F)[i][j] {
                Mx = (*E)[i][j]
                (*P1)[i][j] = (*pE)[i][j]
                (*P2)[i][j] = j
            } else {
                Mx = (*F)[i][j]
                (*P1)[i][j] = i
                (*P2)[i][j] = (*pF)[i][j]
            }

            d1 = (*D)[i][j] + (*D)[i+1][j+1]
            if d1 >= Mx {
                (*D)[i][j] = d1
                (*P1)[i][j] = i+1
                (*P2)[i][j] = j+1
            } else {
                (*D)[i][j] = Mx
            }
        }
    }

    //MxSc := D

    //Bck_Trk
    Dal := make([]IntRow, 2)
    Dal[0] = make(IntRow, 30000)
    Dal[1] = make(IntRow, 30000)

    i := 0
    j := 0
    Alen := 0
    for i < number_of_residues_a && j < number_of_residues_b {
        if i == (*P1)[i][j] && j != (*P2)[i][j] {
            for k := j; k < (*P2)[i][j]; k++ {
                Dal[0][Alen] = -10
                Dal[1][Alen] = k
                Alen += 1
            }
        } else if i != (*P1)[i][j] && j == (*P2)[i][j] {
            for k := i; k < (*P1)[i][j]; k++ {
                Dal[0][Alen] = k
                Dal[1][Alen] = -10
                Alen += 1
            }
        } else if (*P1)[i][j] == i+1 && (*P2)[i][j] == j+1 {
            Dal[0][Alen] = i
            Dal[1][Alen] = j
            Alen += 1
        }
        l := i
        i = (*P1)[i][j]
        j = (*P2)[l][j]
    }

    if i == number_of_residues_a && j < number_of_residues_b {
        for k := j; k < number_of_residues_b; k++ {
            Dal[0][Alen] = -10
            Dal[1][Alen] = k
            Alen += 1
        }
    } else if i < number_of_residues_a && j == number_of_residues_b {
        for k := i; k < number_of_residues_a; k++ {
            Dal[0][Alen] = k
            Dal[1][Alen] = -10
            Alen += 1
        }
    }

    Ial := make([]IntRow, 2)
    for i := 0; i < 2; i++ {
        Ial[i] = make(IntRow, Alen)
        for j := 0; j < Alen; j++ {
            Ial[i][j] = Dal[i][j]
        }
    }

    //Format results
    //NOTE: This has been modified so that iner and map12 now are stored in floats!
    nalign := 0
    sim_tot := 0.0
    iner := make(FloatRow, number_of_residues_a+number_of_residues_b)
    map12 := make(FloatMatrix, 0, number_of_residues_a+number_of_residues_b)
    for j := 0; j < Alen; j++ {
        k1 := Ial[0][j]
        k2 := Ial[1][j]
        iner[j] = 0
        if k1 >= 0 && k2 >= 0 {
            iner[j] = (*equivalence_matrix)[k1][k2]
            if iner[j] > 0.0 {
                nalign += 1
                sim_tot += iner[j]
            }
            //fmt.Printf("REALGN %d %d %d\n", k1, k2, iner[j])
            map12 = append(map12, FloatRow{float64(k1), float64(k2), iner[j]})
        } else if k1 >= 0 {
            //fmt.Printf("REALGN %d - 0\n", k1)
            map12 = append(map12, FloatRow{float64(k1), -1, -1})
        } else if k2 >= 0 {
            //fmt.Printf("REALGN - %d 0\n", k2)
            map12 = append(map12, FloatRow{-1, float64(k2), -1})
        }
    }
    return map12, nalign, sim_tot
}

//FASTA
type FASTASequence struct {
    Label string
    Sequence string
}

func NewScannerLarge(file_path string) (*os.File, *bufio.Scanner) {
    file, err := os.Open(file_path)
    check(err)
    scanner := bufio.NewScanner(file)
    buffer_size := 10*1024*1024 //10 MB buffer
    scanner_buffer := make([]byte, 0, buffer_size) //10 MB buffer
    scanner.Buffer(scanner_buffer, buffer_size)
    return file, scanner
}

func ParseFASTA(path string) []FASTASequence {
    sequences := make([]FASTASequence, 0, 10000)
    label := ""
    buffer := bytes.Buffer{}
    file, scanner := NewScannerLarge(path)
    defer file.Close()
    //Parse sequences delimited by new sequences
    for scanner.Scan() {
        line := strings.TrimSpace(scanner.Text())
        if line == "" {
            continue
        }
        if line[0] == '>' {
            sequence := buffer.String()
            if sequence != "" {
                sequences = append(sequences, FASTASequence{label, sequence})
            }
            buffer.Reset()
            label = line[1:]
        } else {
            buffer.WriteString(line)
        }
    }
    //Parse final sequence
    sequence := buffer.String()
    if sequence != "" {
        sequences = append(sequences, FASTASequence{label, sequence})
    }
    return sequences
}

//Get domains for chain
func get_chain_domains(dash_url string, dash_input DASHInput) []RESTDomain {
    unique_domains := make(map[string]RESTDomain)
    response := http_query("GET",
        fmt.Sprintf(DASHChainDomainsURL(dash_url), dash_input.PDBID), nil)

    defer func() {
        io.Copy(ioutil.Discard, response.Body)
        response.Body.Close()
        response.Close = true
    }()
    scanner := bufio.NewScanner(response.Body)
    for scanner.Scan() {
        json_bytes := scanner.Bytes()
        var domain RESTDomain
        err := json.Unmarshal(json_bytes, &domain)
        check(err)
        if domain.StatusCode != -1 {
            fatal(domain)
        }
        domain.ResidueNumberInts = parse_residue_numbers(domain.ResidueNumbers)
        for _, residue_number := range(domain.ResidueNumberInts) {
            if residue_number >= dash_input.Start && residue_number <= dash_input.End {
                _, exist := unique_domains[domain.DomainID]
                if !exist {
                    unique_domains[domain.DomainID] = domain
                }
            }
        }
    }

    domains := make([]RESTDomain, 0, len(unique_domains))
    for _, domain := range(unique_domains) {
        domains = append(domains, domain)
    }
    return domains
}

//Get chain/domain and self-alignments
func get_chain(dash_url string, pdb_id string) RESTChain {
    var chain RESTChain
    response := http_query("GET", fmt.Sprintf(DASHChainURL(dash_url), pdb_id), nil)
    defer func() {
        io.Copy(ioutil.Discard, response.Body)
        response.Body.Close()
        response.Close = true
    }()
    json_bytes, err := ioutil.ReadAll(response.Body)
    err = json.Unmarshal(json_bytes, &chain)
    check(err)
    if chain.StatusCode != -1 {
        fatal(chain)
    }
    return chain
}

func get_chain_self_alignment(dash_url string, pdb_id string) RESTAlignment {
    chain := get_chain(dash_url, pdb_id)
    var alignment RESTAlignment
    alignment.ID1 = chain.PDBID
    alignment.ID2 = chain.PDBID
    alignment.PRIMS1 = replace_non_standard_residues(chain.Sequence)
    alignment.PRIMS2 = alignment.PRIMS1
    alignment.SECOS1 = strings.Repeat(" ", len(chain.Sequence))
    alignment.SECOS2 = alignment.SECOS1
    alignment.EQUIVALENCE = strings.Repeat("9", len(chain.Sequence))
    return alignment
}

func parse_domain_alignment(json_string string) RESTAlignment {
    var alignment RESTAlignment
    err := json.Unmarshal([]byte(json_string), &alignment)
    check(err)
    if alignment.StatusCode == -1 {
        alignment.PRIMS1 = replace_non_standard_residues(alignment.PRIMS1)
        alignment.PRIMS2 = replace_non_standard_residues(alignment.PRIMS2)
    } else if alignment.StatusCode == 17 {
        alignment.LOWSIMILARITY = true
    } else {
        non_fatal(alignment)
    }
    return alignment
}

func dummy_prims(prims string, start int, end int) (string, int, int) {
    //Create new dummy of prims with all gaps
    new_prims_bytes := make([]byte, len(prims))
    for i, _ := range(new_prims_bytes) {
        new_prims_bytes[i] = '-'
    }
    //Fill in only residues between start and end, save indices to slice later
    start_index := -1
    end_index := -1
    count := 0
    for i := 0; i < len(prims); i++ {
        if count > end {
            break
        }
        if prims[i] != '-' {
            count += 1
        }
        if prims[i] != '-' && count >= start && count <= end {
            new_prims_bytes[i] = prims[i]
        }
        if count == start && start_index == -1 {
            start_index = i
        }
        if count == end && end_index == -1 {
            end_index = i
        }
    }
    return string(new_prims_bytes), start_index, end_index+1
}

func slice_alignment(alignment RESTAlignment,
    query_start int, query_end int, subject_start int, subject_end int) RESTAlignment {
    new_query_prims, start_index, end_index :=
        dummy_prims(alignment.PRIMS1, query_start, query_end)
    new_subject_prims, subject_start_index, subject_end_index :=
        dummy_prims(alignment.PRIMS2, subject_start, subject_end)
    alignment.PRIMS1 = new_query_prims
    alignment.PRIMS2 = new_subject_prims
    if subject_start_index < start_index {
        start_index = subject_start_index
    }
    if subject_end_index > end_index {
        end_index = subject_end_index
    }
    alignment.PRIMS1 = alignment.PRIMS1[start_index:end_index]
    alignment.PRIMS2 = alignment.PRIMS2[start_index:end_index]
    alignment.SECOS1 = alignment.SECOS1[start_index:end_index]
    alignment.SECOS2 = alignment.SECOS2[start_index:end_index]
    alignment.EQUIVALENCE = alignment.EQUIVALENCE[start_index:end_index]
    new_equivalence_bytes := make([]byte, len(alignment.EQUIVALENCE))
    for i := 0; i < len(alignment.PRIMS1); i++ {
        if alignment.PRIMS1[i] != '-' && alignment.PRIMS2[i] != '-' {
            new_equivalence_bytes[i] = alignment.EQUIVALENCE[i]
        } else {
            new_equivalence_bytes[i] = '0'
        }
    }
    alignment.EQUIVALENCE = string(new_equivalence_bytes)
    new_score := 0
    for _, equivalence_byte := range(new_equivalence_bytes) {
        equivalence, err := strconv.Atoi(string(equivalence_byte))
        check(err)
        new_score += equivalence
    }
    alignment.SCORE = new_score
    return alignment
}

//Output formatting functions
func format_alignment_legacy(alignment RESTAlignment) string {
    lines := make([]string, 0, 10)
    if alignment.LOWSIMILARITY {
        lines = append(lines,
            fmt.Sprintf("Query %s Template %s lowsimilarity",
                alignment.ID1, alignment.ID2))
    } else {
        lines = append(lines, fmt.Sprintf("Query %s Template %s",
            alignment.ID1, alignment.ID2))
    }
	lines = append(lines, fmt.Sprintf("QUERY           %s", alignment.PRIMS1))
	lines = append(lines, fmt.Sprintf("QUERY           %s", alignment.SECOS1))
	lines = append(lines, fmt.Sprintf("TEMPL           %s", alignment.PRIMS2))
	lines = append(lines, fmt.Sprintf("TEMPL           %s", alignment.SECOS2))
	lines = append(lines, fmt.Sprintf("Equivalence     %s", alignment.EQUIVALENCE))
    lines = append(lines, "")
    return strings.Join(lines, "\n")
}

func average_over_window(equivalences *[]float64, index int, lookaround int) float64 {
    if lookaround == 0 {
        return (*equivalences)[index]
    }
    start := index - lookaround
    if start < 0 {
        start = 0
    }
    end := index + lookaround + 1
    if end > len(*equivalences) {
        end = len(*equivalences)
    }
    total := 0.0
    count := 0.0
    for i := start; i < end; i++ {
        equivalence := (*equivalences)[i]
        if equivalence > 0 {
            total += equivalence
            count += 1.0
        }
    }
    return total/count
}

func output_alignment_hat3(hat3_file *os.File,
    query_index int, template_index int, query string,
    template string, equivalence_string string, equivalence_threshold float64,
    equivalence_scale float64, minimum_segment_length int,
    equivalence_lookaround int) {

    equivalences := make([]float64, len(query))
    equivalence_mask := make([]bool, len(query))
    //Parse equivalence, enforce threshold, rescale equivalence, multiply by scale
    for i := 0; i < len(query); i++ {
        equivalence_int, err := strconv.Atoi(string(equivalence_string[i]))
        check(err)
        equivalence := float64(equivalence_int)

        if equivalence >= equivalence_threshold {
            equivalence := (equivalence - equivalence_threshold+1) /
                (9 - equivalence_threshold+1) * 9.0
            equivalence = equivalence * equivalence_scale
            equivalences[i] = equivalence
        }
    }

    //Make a mask based on minimum segment length
    last_zero_index := -1
    last_index := len(equivalences) - 1
    for x, equivalence := range(equivalences) {
        //Case where the end of a segment is a 0 equivalence
        if equivalence <= 0 {
            var start_index int
            var length int
            if last_zero_index == -1 {
                //If there is no last zero index
                length = x
                start_index = 0
            } else {
                //If there is
                length = x - last_zero_index - 1
                start_index = last_zero_index + 1
            }
            //Fill mask if long enough
            if length >= minimum_segment_length {
                for start_index < x {
                    equivalence_mask[start_index] = true
                    start_index += 1
                }
            }
            //Set this as the previous index with a zero
            last_zero_index = x
        //Special case for when the alignment ends without a zero equivalence
        } else if x == last_index {
            var start_index int
            var length int
            if last_zero_index == -1 {
                //If there is no last zero index then the length is the whole thing
                length = len(equivalences)
                start_index = 0
            } else {
                //If there is
                length = x - last_zero_index
                start_index = last_zero_index + 1
            }
            //Fill mask if long enough
            if length >= minimum_segment_length {
                for start_index <= x {
                    equivalence_mask[start_index] = true
                    start_index += 1
                }
            }
        }
    }

    //Compute final equivalence and construct line
    lines := make([]string, 0, len(query))
    query_i := -1
    template_i := -1
    for i := 0; i < len(query); i++ {
        if query[i] != '-' {
            query_i += 1
        }
        if template[i] != '-' {
            template_i += 1
        }
        if equivalence_mask[i] {
            windowed_equivalence := average_over_window(&equivalences, i,
                equivalence_lookaround)
            lines = append(lines, fmt.Sprintf("%d %d 1 %0.5f %d %d %d %d k",
                query_index, template_index, windowed_equivalence,
                query_i, query_i, template_i, template_i))
        }
    }

    if len(lines) > 0 {
        fmt.Fprintln(hat3_file, strings.Join(lines, "\n"))
    }
}

func replace_non_standard_residues(sequence string) string {
    sequence = strings.ToUpper(sequence)
    sequence = strings.Replace(sequence, "U", "X", -1)
    sequence = strings.Replace(sequence, "J", "X", -1)
    sequence = strings.Replace(sequence, "O", "X", -1)
    return sequence
}

func filter_sequences_and_hat3(sequence_file_path string, hat3_path string,
    minimum_alignment_percent float64) {
    //Read Sequence Data and Save Indexes
    sequence_data := make([]Sequence, 0, 10000)
    sequence_file, err := os.Open(sequence_file_path)
    check(err)
    scanner := bufio.NewScanner(sequence_file)
    for scanner.Scan() {
        id := strings.TrimSpace(scanner.Text())[1:]
        scanner.Scan()
        sequence := strings.TrimSpace(scanner.Text())
        sequence_data = append(sequence_data, Sequence{id, sequence})
    }
    sequence_file.Close()

    //Read hat3 file to see which ID's were used
    hat3_file, err := os.Open(hat3_path)
    check(err)
    scanner = bufio.NewScanner(hat3_file)
    old_hat3_id_map := make(map[int][]int)
    previous_query_index := -1
    previous_subject_index := -1
    for scanner.Scan() {
        line := strings.TrimSpace(scanner.Text())
        if line == "" {
            continue
        }
        fields := strings.Fields(line)
        query_index, err := strconv.Atoi(fields[0])
        check(err)
        subject_index, err := strconv.Atoi(fields[1])
        check(err)
        if query_index != previous_query_index ||
            subject_index != previous_subject_index {
            old_hat3_id_map[query_index] =
                append(old_hat3_id_map[query_index], subject_index)
            old_hat3_id_map[subject_index] =
                append(old_hat3_id_map[subject_index], query_index)
            previous_query_index = query_index
            previous_subject_index = subject_index
       }
    }
    hat3_file.Close()

    //Append alignment counts to ID's
    for query_index, subject_indexes := range(old_hat3_id_map) {
        sequence_data[query_index].Label += fmt.Sprintf("||%d", len(subject_indexes))
    }

    //Filter least-used sequences from hat3
    alignment_count_cutoff :=
        int(float64(len(old_hat3_id_map))*minimum_alignment_percent/100.0)
    filtered_old_hat3_index_map := make(map[int]bool)
    for query_index, subject_indexes := range(old_hat3_id_map) {
        if len(subject_indexes) >= alignment_count_cutoff {
            filtered_old_hat3_index_map[query_index] = true
            for _, subject_index := range(subject_indexes) {
                if len(old_hat3_id_map[subject_index]) >= alignment_count_cutoff {
                    filtered_old_hat3_index_map[subject_index] = true
                }
            }
        }
    }

    //Make map from new to old indexes
    existing_hat3_indexes := make([]int, 0, len(old_hat3_id_map))
    for i, _ := range(filtered_old_hat3_index_map) {
        existing_hat3_indexes = append(existing_hat3_indexes, i)
    }
    sort.Ints(existing_hat3_indexes)
    hat3_id_map_old_new := make(map[int]int)
    for new_index, old_index := range(existing_hat3_indexes) {
        hat3_id_map_old_new[old_index] = new_index
    }

    //Write new hat3 file
    new_hat3_path := fmt.Sprintf("%s_cleaned", hat3_path)
    new_hat3_file, err := os.Create(new_hat3_path)
    check(err)
    hat3_file, err = os.Open(hat3_path)
    check(err)
    scanner = bufio.NewScanner(hat3_file)
    for scanner.Scan() {
        fields := strings.Fields(strings.TrimSpace(scanner.Text()))
        query_index, err := strconv.Atoi(fields[0])
        check(err)
        subject_index, err := strconv.Atoi(fields[1])
        check(err)
        new_query_index, new_query_index_exist :=
            hat3_id_map_old_new[query_index]
        new_subject_index, new_subject_index_exist :=
            hat3_id_map_old_new[subject_index]
        if new_query_index_exist && new_subject_index_exist {
            fields[0] = fmt.Sprintf("%d", new_query_index)
            fields[1] = fmt.Sprintf("%d", new_subject_index)
            fmt.Fprintln(new_hat3_file, strings.Join(fields, " "))
        }
    }
    hat3_file.Close()
    new_hat3_file.Close()

    //Write Cleaned Sequences
    new_sequence_file, err := os.Create(sequence_file_path)
    check(err)
    defer new_sequence_file.Close()
    for _, old_index := range(existing_hat3_indexes) {
        sequence := sequence_data[old_index]
        sequence_lines := fmt.Sprintf(">%s\n%s", sequence.Label, sequence.Sequence)
        fmt.Fprintln(new_sequence_file, sequence_lines)
    }
    //Move new hat3 to old hat3
    err = os.Rename(new_hat3_path, hat3_path)
    check(err)
}

func main() {
	fmt.Println("------------------")
	fmt.Println("MAFFT-DASH Client v", VersionNumber)
	fmt.Println("------------------")
    //Parse flags
    var help bool
    var input_path string
    var hat3_output_path string
    var alignments_output_path string
    var sequences_output_path string
    var slice bool
    var filter float64
    var equivalence_threshold float64
    var equivalence_scale float64
    var blosum_alpha float64
    var structure_only bool
    var minimum_segment_length int
    var equivalence_lookaround int
    var template_list_path string
    var dash_url string
    flag.BoolVar(&help, "help", false, "Display this help message.")
    flag.StringVar(&input_path, "i", "",
        fmt.Sprintf("Path to FASTA sequence file. (REQUIRED!) (Limit of %d sequences)",
            SequenceLimit))
    flag.StringVar(&dash_url, "url", "https://sysimm.org/dash/REST1.0/",
        "URL for DASH REST service.")
    flag.StringVar(&template_list_path, "templates", "",
        "Path to explicit list of DASH templates (Debug-only).")
    flag.BoolVar(&slice, "slice", true,
        "Slice alignments/sequences according to start/end positions.")
    flag.Float64Var(&filter, "filter", 22.5,
        "Filter sequences/hat3 where sequence is aligned to less than this percentage of the inputs.")
    flag.Float64Var(&equivalence_threshold, "threshold", 1.0,
        "Only use equivalence values >= this value for hat3.")
    flag.Float64Var(&equivalence_scale, "scale", 1.0,
        "Multiply equivalence values by this value before outputting to hat3.")
    flag.IntVar(&minimum_segment_length, "length", 5,
        "Only use equivalences when the number of consecutive residues >= the threshold are >= this length.")
    flag.IntVar(&equivalence_lookaround, "lookaround", 0,
        "Output the average of this many surrounding values instead of raw equivalence. (default 0)")
    flag.StringVar(&hat3_output_path, "hat3", "./hat3", "Output path for hat3 file.")
    flag.StringVar(&alignments_output_path, "alignments", "./dash_alignments",
        "Output path for raw alignments for debug purposes.")
    flag.StringVar(&sequences_output_path, "sequences", "./dash_sequences.fa",
        "Output path for template sequences for debug purposes.")
    flag.Float64Var(&blosum_alpha, "alpha", 0.75, "Background sequence BLOSUM multiplier.")
    flag.BoolVar(&structure_only, "structure-only", false,
        "Output alignments with only residues which exist in the structure. Sets BLOSUM Alpha to 0.0")
    flag.Parse()
    if input_path == "" && template_list_path == "" {
        fmt.Println(
            "Please submit a FASTA sequence file with -i or a template list with -templates.")
        fmt.Println("------------------")
        flag.PrintDefaults()
        os.Exit(0)
    }
    if help {
        flag.PrintDefaults()
        os.Exit(0)
    }
    if structure_only {
        blosum_alpha = 0.0
    }
    fmt.Println("Querying from", dash_url)

    //Parse sequences and query DASH for representatives
    sequences := []FASTASequence{}
    if input_path != "" {
        sequences = ParseFASTA(input_path)
        if len(sequences) > SequenceLimit {
            fatal(fmt.Sprintf("Number of sequences greater than sequence limit of %d.",
                SequenceLimit))
        }
    }
    template_list_buffer := bytes.Buffer{}
    if template_list_path != "" {
        //Parse template list into JSON so that it's the same format as what the
        //server would return
        fmt.Println("Using local template list from", template_list_path)
        template_list, err := os.Open(template_list_path)
        check(err)
        scanner := bufio.NewScanner(template_list)
        i := 0
        for scanner.Scan() {
            i += 1
            line := strings.TrimSpace(scanner.Text())
            if line == "" {
                continue
            }
            split_line := strings.Split(line, "||")
            if len(split_line) != 3 {
                fatal(fmt.Sprintf("Line %d contains wrong number of fields.", i))
            }
            id := split_line[0]
            start, err := strconv.Atoi(split_line[1])
            check(err)
            end, err := strconv.Atoi(split_line[2])
            check(err)
            object := RESTSearch{-1, "", id, start, end}
            object_bytes, err := json.Marshal(object)
            check(err)
            _, err = template_list_buffer.Write(object_bytes)
            check(err)
            _, err = template_list_buffer.WriteString("\n")
            check(err)
        }
        template_list.Close()
    } else if input_path != "" {
        //Write request body for template selection
        fmt.Println("Building query for template selection.")
        request_body := bytes.Buffer{}
        for i, sequence := range(sequences) {
            //Use index if labels don't exist
            label := fmt.Sprintf("%d", i)
            if sequence.Label != "" {
                label = sequence.Label
            }
            _, err := request_body.WriteString(fmt.Sprintf(">%s\n%s\n",
                label, sequence.Sequence))
            check(err)
        }
        //Send query
        fmt.Println("Sending query for template selection...")
        response := http_query("POST", DASHChainSearchURL(dash_url), &request_body)
        _, err := io.Copy(&template_list_buffer, response.Body)
        check(err)
        io.Copy(ioutil.Discard, response.Body)
        response.Body.Close()
        fmt.Println("Waiting for response from server...")
    }

    //Parse results for template selection
    dash_inputs := make([]DASHInput, 0, 10*len(sequences))
    id_map := make(map[string]bool)
    scanner := bufio.NewScanner(&template_list_buffer)
    for scanner.Scan() {
        json_bytes := scanner.Bytes()
        var result RESTSearch
        err := json.Unmarshal(json_bytes, &result)
        check(err)
        if result.StatusCode != -1 {
            fatal(result)
        }
        var dash_input DASHInput
        dash_input.FullID =
            fmt.Sprintf("%s||%d||%d", result.ID, result.Start, result.End)
        dash_input.PDBID = result.ID
        dash_input.Start = result.Start
        dash_input.End = result.End
        _, exists := id_map[dash_input.FullID]
        if !exists {
            dash_inputs = append(dash_inputs, dash_input)
            id_map[dash_input.FullID] = true
        }
    }

    sort.Slice(dash_inputs, func(i, j int) bool {
        return dash_inputs[i].FullID < dash_inputs[j].FullID
    })

    //Split chains into domains
    dash_input_map := make(map[string][]DASHInput)
    dash_domain_id_map := make([]string, 0, 10000)
    for i, dash_input := range(dash_inputs) {
        if i % 25 == 0 {
            percent := i*100/len(dash_inputs)
            fmt.Printf("Querying DASH for domains - [%d%%]\n", percent)
        }
        dash_input.Hat3Index = i
        dash_input.Domains = get_chain_domains(dash_url, dash_input)
        for _, domain := range(dash_input.Domains) {
            _, exist := dash_input_map[domain.DomainID]
            if !exist {
                dash_domain_id_map = append(dash_domain_id_map, domain.DomainID)
            }
            dash_input_map[domain.DomainID] =
                append(dash_input_map[domain.DomainID], dash_input)
        }
        dash_inputs[i] = dash_input
    }
    fmt.Printf("Querying DASH for domains - [100%%]\n")

    //Open output files
    sequences_output_file, err := os.Create(sequences_output_path)
    check(err)
    alignments_output_file, err := os.Create(alignments_output_path)
    check(err)
    hat3_output_file, err := os.Create(hat3_output_path)
    check(err)

    //Get sequences and self-alignments
    for i, dash_input := range(dash_inputs) {
        if i % 25 == 0 {
            percent := i*100/len(dash_inputs)
            fmt.Printf("Querying DASH for sequences - [%d%%]\n", percent)
        }
        raw_alignment := get_chain_self_alignment(dash_url, dash_input.PDBID)
        alignment := raw_alignment
        if slice {
            alignment = slice_alignment(alignment,
                dash_input.Start, dash_input.End,
                dash_input.Start, dash_input.End)
        }
        dash_inputs[i].Sequence = alignment.PRIMS1
        fmt.Fprintln(sequences_output_file,
            fmt.Sprintf(">DASH_%s\n%s",
            dash_input.FullID, alignment.PRIMS1))
        alignment.ID1 = dash_input.FullID
        alignment.ID2 = dash_input.FullID
        fmt.Fprintln(alignments_output_file,
            format_alignment_legacy(alignment))
    }
    fmt.Printf("Querying DASH for sequences - [100%%]\n")

    //Batch alignment downloads in groups based on alignment limit
    alignment_map := make(map[string]RESTAlignment)
    request_body := bytes.Buffer{}
    alignment_count := 0
    chunk_index := 0
    number_of_domains := len(dash_domain_id_map)
    number_of_domain_alignments := (number_of_domains*(number_of_domains-1)/2)
    for x := 0; x < len(dash_domain_id_map); x++ {
        for y := x + 1; y < len(dash_domain_id_map); y++ {
            _, err := request_body.WriteString(fmt.Sprintf("%s_%s\n",
                dash_domain_id_map[x], dash_domain_id_map[y]))
            check(err)
            alignment_count += 1
            if alignment_count >= AlignmentLimit {
                current_progress := chunk_index*AlignmentLimit
                current_percent := current_progress*100/number_of_domain_alignments
                fmt.Printf("Downloading alignments - [%d%%]\n", current_percent)
                //Submit to server
                response := http_query("POST",
                    DASHDomainAlignmentURL(dash_url), &request_body)
                scanner = bufio.NewScanner(response.Body)
                //Parse results
                for scanner.Scan() {
                    alignment := parse_domain_alignment(scanner.Text())
                    if alignment.StatusCode == -1 {
                        alignment_id := alignment.ID1 + "_" + alignment.ID2
                        alignment_map[alignment_id] = alignment
                    }
                }
                //Reset
                chunk_index += 1
                io.Copy(ioutil.Discard, response.Body)
                response.Body.Close()
                alignment_count = 0
                request_body.Reset()
            }
        }
    }
    //Get remaining alignments
    if alignment_count > 0 {
        current_progress := chunk_index*AlignmentLimit
        current_percent := current_progress*100/number_of_domain_alignments
        fmt.Printf("Downloading alignments - [%d%%]\n", current_percent)
        //Submit to server
        response := http_query("POST",
            DASHDomainAlignmentURL(dash_url), &request_body)
        scanner = bufio.NewScanner(response.Body)
        //Parse results
        for scanner.Scan() {
            alignment := parse_domain_alignment(scanner.Text())
            if alignment.StatusCode == -1 {
                alignment_id := alignment.ID1 + "_" + alignment.ID2
                alignment_map[alignment_id] = alignment
            }
        }
        io.Copy(ioutil.Discard, response.Body)
        response.Body.Close()
    }
    fmt.Printf("Downloading alignments - [100%%]\n")

    //Combine domain alignments into full chain alignments
    number_of_chains := len(dash_inputs)
    number_of_chain_alignments := (number_of_chains*(number_of_chains-1)/2)
    alignment_count = 0
    for x := 0; x < len(dash_inputs); x++ {
        for y := x + 1; y < len(dash_inputs); y++ {
            if alignment_count % AlignmentLimit == 0 {
                current_percent := alignment_count*100/number_of_chain_alignments
                fmt.Printf("Combining domain alignments - [%d%%]\n", current_percent)
            }
            alignment_count += 1
            dash_input_a := dash_inputs[x]
            dash_input_b := dash_inputs[y]
            size_a := len(dash_input_a.Sequence)
            size_b := len(dash_input_b.Sequence)
            //Initialize equivalence matrix with BLOSUM
            equivalence_matrix := InitializeFloatMatrix(size_a, size_b)
            low_similarity := true
            for x, query_residue := range dash_input_a.Sequence {
                for y, template_residue := range dash_input_b.Sequence {
                    blosum_score :=
                        float64(BLOSUM62[byte(query_residue)][byte(template_residue)])
                    if blosum_score > 0 {
                        //Re-scale BLOSUM scores to RASH's 0-9.99 scale
                        blosum_score = blosum_score / float64(BLOSUM62Max) * 9.999
                        //Multiply by a factor to control the influence
                        (*equivalence_matrix)[x][y] = blosum_score * blosum_alpha
                    }
                }
            }
            for _, domain_a := range(dash_input_a.Domains) {
                for _, domain_b := range(dash_input_b.Domains) {
                    //Fill matrix with equivalence data for realignment
                    alignment_id := domain_a.DomainID + "_" + domain_b.DomainID
                    alignment, exist := alignment_map[alignment_id]
                    if !exist {
                        alignment_id = domain_b.DomainID + "_" + domain_a.DomainID
                        alignment, exist = alignment_map[alignment_id]
                        if exist {
                            alignment.Reverse()
                        }
                    }

                    if exist {
                        low_similarity = false
                        index_a := -1
                        index_b := -1
                        for i := 0; i < len(alignment.PRIMS1); i++ {
                            if alignment.PRIMS1[i] != '-' { index_a += 1 }
                            if alignment.PRIMS2[i] != '-' { index_b += 1 }
                            if alignment.PRIMS1[i] != '-' && alignment.PRIMS2[i] != '-' {
                                matrix_index_a :=
                                    domain_a.ResidueNumberInts[index_a] - dash_input_a.Start
                                matrix_index_b :=
                                    domain_b.ResidueNumberInts[index_b] - dash_input_b.Start
                                equivalence, err :=
                                    strconv.ParseFloat(string(alignment.EQUIVALENCE[i]), 64)
                                check(err)
                                if matrix_index_a >= 0 && matrix_index_a < size_a &&
                                    matrix_index_b >= 0 && matrix_index_b < size_b {
                                    (*equivalence_matrix)[matrix_index_a][matrix_index_b] =
                                        equivalence
                                }
                            }
                         }
                    }
                }
            }

            //Construct alignment from realignment
            aligned_matrix, _, _ := AlignMatrix(equivalence_matrix)
            PRIMS1_bytes := []byte(strings.Repeat("-", len(aligned_matrix)))
            PRIMS2_bytes := []byte(strings.Repeat("-", len(aligned_matrix)))
            EQUIVALENCE_bytes := []byte(strings.Repeat("0", len(aligned_matrix)))
            for i, row := range(aligned_matrix) {
                index_a := int(row[0])
                index_b := int(row[1])
                equivalence := strconv.Itoa(int(row[2]))
                if index_a != -1 && index_b != -1 {
                    EQUIVALENCE_bytes[i] = equivalence[0]
                }
                if index_a != -1 {
                    PRIMS1_bytes[i] = dash_input_a.Sequence[index_a]
                }
                if index_b != -1 {
                    PRIMS2_bytes[i] = dash_input_b.Sequence[index_b]
                }
            }
            var alignment RESTAlignment
            alignment.ID1 = dash_input_a.FullID
            alignment.ID2 = dash_input_b.FullID
            alignment.PRIMS1 = string(PRIMS1_bytes)
            alignment.PRIMS2 = string(PRIMS2_bytes)
            alignment.SECOS1 = strings.Repeat(" ", len(aligned_matrix))
            alignment.SECOS2 = alignment.SECOS1
            alignment.EQUIVALENCE = string(EQUIVALENCE_bytes)
            alignment.LOWSIMILARITY = low_similarity

            //Output alignment and hat3
            fmt.Fprintln(alignments_output_file,
                format_alignment_legacy(alignment))
            if !alignment.LOWSIMILARITY {
                output_alignment_hat3(hat3_output_file,
                    dash_input_a.Hat3Index, dash_input_b.Hat3Index,
                    alignment.PRIMS1, alignment.PRIMS2,
                    alignment.EQUIVALENCE, equivalence_threshold,
                    equivalence_scale, minimum_segment_length,
                    equivalence_lookaround)
            }
        }
    }
    fmt.Printf("Combining domain alignments - [100%%]\n")

    //Filter unused sequences from hat3 and sequences file which are below the threshold
    fmt.Println("Filtering structural restraint(hat3) file...")
    sequences_output_file.Close()
    hat3_output_file.Close()
    filter_sequences_and_hat3(sequences_output_path, hat3_output_path, filter)

    //Combine original sequences with DASH sequences
    fmt.Println("Combining original sequences with DASH sequences...")
    dash_sequences := ParseFASTA(sequences_output_path)
    final_sequence_file, err := os.Create(sequences_output_path)
    check(err)
    for _, sequence := range(dash_sequences) {
        fmt.Fprintln(final_sequence_file,
            fmt.Sprintf(">%s\n%s", sequence.Label, sequence.Sequence))
    }
    for _, sequence := range(sequences) {
        fmt.Fprintln(final_sequence_file,
            fmt.Sprintf(">%s\n%s", sequence.Label, sequence.Sequence))
    }
    final_sequence_file.Close()

    //Final user output
	fmt.Println("------------------")
    fmt.Println("Ready to run MAFFT:")
    fmt.Println("  mafft --seedtable", hat3_output_path,
        "--localpair", "--maxiterate 100", sequences_output_path)
}
